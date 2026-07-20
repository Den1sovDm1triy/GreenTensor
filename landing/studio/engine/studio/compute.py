# SPDX-License-Identifier: MIT

"""compute — расчётный движок GreenTensor Studio (поверх ПУБЛИЧНОГО API библиотеки).

Studio — отдельный потребитель библиотеки. ОБЛАСТЬ — только точные аналитические
ТФГ-решения: слоистая СФЕРА (каноническое ядро ``01_sphere.py`` через
``green_tensor.SphereSolver``, фасад ``sphere_core``), бесконечный слоистый ЦИЛИНДР
(``LayeredCylinderSolver``, точная 2D-аналитика) и КЛАСТЕР невзаимопересекающихся
сфер (``Cluster``/GMM, теорема сложения Крузана–Стейна). Ближнепольные тепловые
карты — через ``Cluster.scattered_field`` (для сфер и их кластеров).

Вход — сцена (см. presets.py), выход — словарь серий для графиков (JSON).
"""
from __future__ import annotations

import cmath
import math

import numpy as np
import scipy.special as _sp

from .library_bridge import load as _load_library

DB_FLOOR = -60.0
MAX_TOCH = 200
MAX_SWEEP_POINTS = 400
MAX_HEATMAP = 64          # точек на сторону near-field карты (по умолчанию)
HEATMAP_RES_MIN = 32      # минимум сетки карты
HEATMAP_RES_MAX = 200     # жёсткий предел сетки (защита производительности)
N_THETA = 720             # узлов диаграммы (полный круг)
FAR_R_FACTOR = 200.0      # радиус выборки дальнего поля = FAR_R_FACTOR / k


class ComputeError(ValueError):
    """Ошибка валидации входной сцены (текст показывается в UI)."""


# --------------------------------------------------------------------------- #
# Утилиты
# --------------------------------------------------------------------------- #
def _finite(arr):
    return [float(v) if np.isfinite(v) else None for v in np.asarray(arr, float).ravel()]


def _db(values, ref=None):
    """dB-нормировка к общему максимуму ref (или к собственному максимуму)."""
    v = np.abs(np.asarray(values, float))
    vmax = float(ref) if ref is not None else (float(v[np.isfinite(v)].max()) if v.size else 0.0)
    if vmax <= 0:
        return np.full_like(v, DB_FLOOR)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = 20.0 * np.log10(v / vmax)
    return np.maximum(d, DB_FLOOR)


def _cnum(re, im=0.0):
    return complex(float(re), float(im))


def _layers(body):
    """eps/mu/a_norm (изнутри наружу) из слоёв тела сцены."""
    raw = body.get("layers") or []
    if not raw:
        raise ComputeError("У тела не заданы слои/материал.")
    a_norm, eps, mu = [], [], []
    for i, L in enumerate(raw, 1):
        try:
            a_norm.append(float(L["a_rel"]))
            eps.append(_cnum(L.get("eps_re", 1.0), L.get("eps_im", 0.0)))
            mu.append(_cnum(L.get("mu_re", 1.0), L.get("mu_im", 0.0)))
        except (KeyError, TypeError, ValueError):
            raise ComputeError(f"Слой {i}: некорректные значения.")
    if any(b <= a for a, b in zip(a_norm, a_norm[1:])):
        raise ComputeError("Границы слоёв a_rel должны строго возрастать (внешняя = 1).")
    return a_norm, eps, mu


def _radiation(scene):
    r = scene.get("radiation") or {}
    try:
        k = float(r.get("k", 0.0))
    except (TypeError, ValueError):
        raise ComputeError("k: некорректное значение.")
    if not (k > 0 and math.isfinite(k)):
        raise ComputeError("k должно быть > 0.")
    pol_mode = r.get("polarization", "linear")
    if pol_mode not in ("linear", "circular"):
        raise ComputeError("polarization: linear|circular.")
    toch = r.get("toch")
    if toch is not None:
        toch = int(toch)
        if not (1 <= toch <= MAX_TOCH):
            raise ComputeError(f"toch: допустимо 1..{MAX_TOCH}.")
    khat = np.asarray(r.get("khat", [0.0, 0.0, 1.0]), float)
    khat = khat / (np.linalg.norm(khat) or 1.0)
    pol = np.asarray(r.get("pol", [1.0, 0.0, 0.0]), float)
    return {
        "k": k, "polarization": pol_mode, "problem": r.get("problem", "diffraction"),
        "phi_deg": float(r.get("phi_deg", 0.0)), "toch": toch,
        "khat": khat, "pol": pol,
    }


# --------------------------------------------------------------------------- #
# Построение рассеивателей (для кластера/GMM и ближнего поля)
# --------------------------------------------------------------------------- #
def _build_scatterer(body, gt):
    """Тело сцены -> рассеиватель публичного API (для Cluster / near-field).

    Поддерживаются только СФЕРЫ (точное ядро 01_sphere); кластер собирается строгой
    теоремой сложения GMM. Бесконечный цилиндр считается отдельно как одиночное тело.
    """
    btype = body.get("type", "sphere")
    pos = body.get("position", [0.0, 0.0, 0.0])
    a_norm, eps, mu = _layers(body)

    if btype == "sphere":
        radius = float(body.get("radius", 1.0))
        return gt.SphereSolver(radius, eps, a_norm=a_norm,
                               miy=mu, position=pos).as_scatterer()
    if btype == "cylinder_inf":
        raise ComputeError("Бесконечный цилиндр — 2D-задача; поддерживается только как "
                           "одиночное тело, не в кластере сфер.")
    raise ComputeError(
        f"Тип тела {btype!r} не поддерживается. Точная ТФГ-аналитика доступна для сферы "
        f"(в т.ч. в кластере сфер) и одиночного бесконечного цилиндра."
    )


def _nmax_for(bodies, k):
    """Оценка nmax для GMM по РЕАЛЬНЫМ радиусам сфер (минимум nmax=4, максимум 16)."""
    rmax = max((float(b.get("radius", 0.0)) for b in bodies), default=0.0)
    return int(min(max(math.ceil(k * rmax + 4.0 * (k * rmax) ** (1 / 3) + 2.0), 4), 16))


# --------------------------------------------------------------------------- #
# Полусфера на проводящем экране (метод зеркальных изображений)
# --------------------------------------------------------------------------- #
HEMI_N_THETA = 720          # узлов по кругу; сетка сдвинута на полшага (углы ±90° не на узлах)


MAX_HEMI_FEEDS = 32


def _hemisphere_cfg(body):
    """Конфиг полусферы тела или None.

    Возвращает {'feeds': [{'offset_deg', 'amp', 'phase_deg'}, ...]}.
    Одиночный облучатель задаётся либо feeds из одного элемента, либо
    (обратная совместимость) полем feed_offset_deg. Линейка облучателей с
    амплитудными весами формирует профилированные ДН (например, csc²).
    """
    h = body.get("hemisphere") or {}
    if not h.get("enabled"):
        return None
    raw = h.get("feeds")
    if not raw:
        raw = [{"offset_deg": h.get("feed_offset_deg", 0.0), "amp": 1.0, "phase_deg": 0.0}]
    if len(raw) > MAX_HEMI_FEEDS:
        raise ComputeError(f"Полусфера: не более {MAX_HEMI_FEEDS} облучателей.")
    feeds = []
    for i, f in enumerate(raw, 1):
        try:
            off = float(f.get("offset_deg", 0.0))
            amp = float(f.get("amp", 1.0))
            ph = float(f.get("phase_deg", 0.0))
        except (TypeError, ValueError, AttributeError):
            raise ComputeError(f"Полусфера: облучатель {i} — некорректные значения.")
        if not (math.isfinite(off) and abs(off) <= 80.0):
            raise ComputeError(f"Полусфера: облучатель {i} — смещение θ′ в пределах ±80°.")
        if not (math.isfinite(amp) and 0.0 <= amp <= 100.0):
            raise ComputeError(f"Полусфера: облучатель {i} — амплитуда в пределах 0..100.")
        if not math.isfinite(ph):
            raise ComputeError(f"Полусфера: облучатель {i} — некорректная фаза.")
        feeds.append({"offset_deg": off, "amp": amp, "phase_deg": ph})
    if not any(f["amp"] > 0 for f in feeds):
        raise ComputeError("Полусфера: хотя бы один облучатель должен иметь амплитуду > 0.")
    return {"feeds": feeds}


def _hemisphere_pattern(body, rad, gt):
    """ДН полусферической линзы на бесконечном PEC-экране (image theory).

    Полусфера на экране эквивалентна полной сфере + зеркальному образу источника:
        E_θ(θ) = f_θ(θ−θ′) + f_θ(π−θ−θ′),   E_φ(θ) = f_φ(θ−θ′) − f_φ(π−θ−θ′),
    где θ′ — смещение облучателя от нормали к экрану, θ — угол наблюдения от
    нормали (+z), f — комплексная диаграмма полной сферы. Знаки закреплены
    граничным условием E_tan=0 на экране (E_φ(±90°)≡0 при любом θ′). Комбинация
    эквивалентна отбору сферических гармоник по чётности: τ_n(π−θ)=(−1)ⁿτ_n(θ),
    π_n(π−θ)=(−1)ⁿ⁺¹π_n(θ), поэтому TE- и TM-семейства входят с разной чётностью.
    Смещённый облучатель моделируется наклонным падением плоской волны
    (взаимность); поворот диаграммы точен благодаря сферической симметрии линзы.
    Поле определено в верхнем полупространстве |θ|≤90°; ниже экрана тождественно 0.
    Энергобаланс проверен: мощность зеркального члена в верхней полусфере в
    точности равна мощности прямого поля, уходившей в нижнюю (P_i_up == P_d_down).
    """
    if rad["polarization"] != "linear":
        raise ComputeError("Полусфера на PEC-экране: поддерживается линейная поляризация.")
    radius = float(body.get("radius", 1.0))
    a_norm, eps, mu = _layers(body)
    feeds = _hemisphere_cfg(body)["feeds"]
    k = rad["k"]
    mie = gt.SphereSolver(radius, eps, a_norm=a_norm, miy=mu).mie(k, toch=rad["toch"])
    Mn, Nn = mie.coefficients()
    toch = len(Mn)
    ns = np.arange(1, toch + 1)
    w = ((2 * ns + 1) / (ns * (ns + 1))) * ((-1.0) ** ns)

    def base(u):
        """Комплексные f_θ(u), f_φ(u) полной сферы; u — произвольный угол (рад).

        Чётное 2π-периодическое продолжение среза диаграммы через
        w = arccos(cos u): при переходе через полюс знак присоединённого полинома
        Лежандра компенсируется переворотом базисных векторов среза, поэтому
        компоненты продолжаются чётно и «инверсия Лежандра» не нужна.
        """
        ww = np.arccos(np.clip(np.cos(np.asarray(u, float)), -1.0, 1.0))
        ww = np.clip(ww, 1e-4, math.pi - 1e-4)
        ct, st = np.cos(ww), np.sin(ww)
        ft = np.zeros(ww.shape, complex)
        fp = np.zeros(ww.shape, complex)
        for i in range(toch):
            n = i + 1
            L0 = _sp.lpmv(0, n, ct)
            L1 = _sp.lpmv(1, n, ct)
            L2 = _sp.lpmv(2, n, ct) if n >= 2 else 0.0
            pii = L1 / st
            tay = 0.5 * (L2 - n * (n + 1) * L0)
            ft = ft + w[i] * (tay * Mn[i] - pii * Nn[i])
            fp = fp + w[i] * (pii * Mn[i] - tay * Nn[i])
        return ft, fp

    theta_deg = (np.arange(HEMI_N_THETA) + 0.5) * (360.0 / HEMI_N_THETA)
    th = np.radians(theta_deg)
    # когерентная сумма по линейке облучателей: каждый — прямой луч + зеркальный образ
    Et = np.zeros(th.shape, complex)
    Ep = np.zeros(th.shape, complex)
    for f in feeds:
        if f["amp"] <= 0.0:
            continue
        wgt = f["amp"] * cmath.exp(1j * math.radians(f["phase_deg"]))
        tp = math.radians(f["offset_deg"])
        ft_d, fp_d = base(th - tp)                # прямой (смещённый) облучатель
        ft_i, fp_i = base(math.pi - th - tp)      # зеркальный образ под экраном
        Et = Et + wgt * (ft_d + ft_i)             # нормальная компонента: образ с «+»
        Ep = Ep + wgt * (fp_d - fp_i)             # касательная: образ с «−» (ноль на экране)
    upper = np.cos(th) >= 0.0                     # верхнее полупространство |θ| ≤ 90°
    Et = np.where(upper, Et, 0.0)
    Ep = np.where(upper, Ep, 0.0)
    aEt, aEp = np.abs(Et), np.abs(Ep)
    ref = max(float(aEt.max()), float(aEp.max())) or 1.0
    series = [
        {"name": "E_θ (верт. поляризация)", "E": _finite(aEt), "dB": _finite(_db(aEt, ref))},
        {"name": "E_φ (гориз. поляризация)", "E": _finite(aEp), "dB": _finite(_db(aEp, ref))},
    ]
    return {"theta_deg": _finite(theta_deg), "series": series, "polar": True,
            "plane": "E/H", "hemisphere": {"feeds": feeds},
            "note": "θ отсчитывается от нормали к экрану; ниже экрана (|θ|>90°) поле = 0"}


# --------------------------------------------------------------------------- #
# Диаграмма направленности (ДН)
# --------------------------------------------------------------------------- #
def _sphere_pattern(body, rad, gt):
    """ДН одиночной сферы каноном 01_sphere.py (через SphereSolver)."""
    radius = float(body.get("radius", 1.0))
    a_norm, eps, mu = _layers(body)
    solver = gt.SphereSolver(radius, eps, a_norm=a_norm, miy=mu)
    k = rad["k"]
    # стандартный угол рассеяния θ_std (0=вперёд, 180=обратное), как у кластера;
    # 01_sphere использует обратную конвенцию → запрашиваем θ_01 = π − θ_std
    theta_std = np.linspace(0.0, 2.0 * math.pi, N_THETA, endpoint=False) + 1e-4
    theta = (math.pi - theta_std) % (2.0 * math.pi)
    series = []
    prob = rad["problem"]                                 # diffraction | antenna (источник на сфере)
    if rad["polarization"] == "circular":
        pat = solver.pattern(k, theta, polarization="circular", problem=prob)
        ref = max(np.max(np.abs(pat["E_op"])), np.max(np.abs(pat["E_kp"])))
        for key, name in (("E_op", "Осн. поляризация"), ("E_kp", "Кросс-поляризация")):
            series.append({"name": name, "E": _finite(pat[key]), "dB": _finite(_db(pat[key], ref))})
    else:
        pe = solver.pattern(k, theta, polarization="linear", phi=0.0, problem=prob)["E_theta"]
        ph = solver.pattern(k, theta, polarization="linear", phi=math.pi / 2, problem=prob)["E_theta"]
        ref = max(np.max(np.abs(pe)), np.max(np.abs(ph)))
        series.append({"name": "E-плоскость (φ=0°)", "E": _finite(pe), "dB": _finite(_db(pe, ref))})
        series.append({"name": "H-плоскость (φ=90°)", "E": _finite(ph), "dB": _finite(_db(ph, ref))})
        phi = math.radians(rad["phi_deg"])
        if abs(rad["phi_deg"]) > 1e-6 and abs(rad["phi_deg"] - 90.0) > 1e-6:
            pc = solver.pattern(k, theta, polarization="linear", phi=phi, problem=prob)["E_theta"]
            series.append({"name": f"φ={rad['phi_deg']:.0f}°",
                           "E": _finite(pc), "dB": _finite(_db(pc, ref))})
    return {"theta_deg": _finite(np.degrees(theta_std)), "series": series,
            "polar": True, "plane": "E/H"}


def _cluster_pattern(scatterers, rad, gt, nmax, cf=None):
    """ДН произвольного кластера: |E| дальнего поля в E- и H-плоскостях. Угол θ отсчитывается
    от НАПРАВЛЕНИЯ ПАДЕНИЯ (θ=0 — рассеяние ВПЕРЁД, θ=180 — обратное), как у сферы — иначе для
    непродольного падения метки «вперёд/назад» были бы неверны. E-плоскость содержит вектор
    поляризации ê, H-плоскость — k̂×ê. cf — предрешённые коэф. (если переданы — без пересолва)."""
    k = rad["k"]
    theta = np.linspace(0.0, 2.0 * math.pi, 360, endpoint=False) + 1e-4
    R = FAR_R_FACTOR / k
    cl = gt.Cluster(scatterers)
    if cf is None:
        _, cf, _ = cl.solve(k, tuple(rad["khat"]), tuple(rad["pol"]), nmax)

    khat = np.asarray(rad["khat"], float)
    khat = khat / (np.linalg.norm(khat) or 1.0)          # вперёд = +k̂
    pol = np.asarray(rad["pol"], float)
    pe = pol - (pol @ khat) * khat                        # поляризация ⊥ k̂ (E-плоскость)
    if np.linalg.norm(pe) < 1e-9:                         # вырожденно: ê ∥ k̂ → любой перпендикуляр
        pe = np.cross(khat, [1.0, 0, 0] if abs(khat[0]) < 0.9 else [0, 0, 1.0])
    pe = pe / np.linalg.norm(pe)
    ph = np.cross(khat, pe)                               # H-плоскость, |ph|=1

    def plane(axis):                                      # n̂(θ)=cosθ·k̂ + sinθ·axis; θ=0 — вперёд
        n = np.outer(np.cos(theta), khat) + np.outer(np.sin(theta), axis)
        E = cl.scattered_field(cf, k, R * n)
        return np.linalg.norm(E, axis=1) * R             # |f(θ)| ~ |E|·r

    e_plane = plane(pe)
    h_plane = plane(ph)
    ref = max(np.max(e_plane), np.max(h_plane))
    series = [
        {"name": "E-плоскость (xz)", "E": _finite(e_plane), "dB": _finite(_db(e_plane, ref))},
        {"name": "H-плоскость (yz)", "E": _finite(h_plane), "dB": _finite(_db(h_plane, ref))},
    ]
    return {"theta_deg": _finite(np.degrees(theta)), "series": series,
            "polar": True, "plane": "E/H"}


# --------------------------------------------------------------------------- #
# ЭПР: бистатическая σ(θ) при фикс. падении и моностатическая σ(ракурс).
# Нормировка σ/λ² (λ=2π/k0) — безразмерная и однозначная для любой геометрии.
# --------------------------------------------------------------------------- #
MONO_STEP_DEG = 2.5            # шаг развёртки ракурса в моностатической ЭПР


def _sigma_db(sig):
    return 10.0 * np.log10(np.maximum(np.asarray(sig, dtype=float), 1e-300))


def _bistatic_from_pattern(pattern, const):
    """σ/λ² по плоскостям из амплитуд ДН: σ/λ² = const·|E|²."""
    series = []
    for s in pattern["series"]:
        E = np.array([v if v is not None else 0.0 for v in s["E"]], dtype=float)
        sig = const * E * E
        series.append({"name": s["name"], "sigma": _finite(sig), "dB": _finite(_sigma_db(sig))})
    return {"angle_deg": list(pattern["theta_deg"]), "series": series,
            "unit": "σ/λ²", "available": True}


def _bistatic_cluster(pattern, k):
    # σ = 4π|f|², f = |E|·R (амплитуда ДН); σ/λ² = |f|²·k²/π
    return _bistatic_from_pattern(pattern, const=k * k / math.pi)


def _bistatic_sphere(pattern, k, radius, q_back):
    # якорь по обратному рассеянию: σ_back/λ² = q_back·πR²/λ² = q_back·(kR)²/(4π)
    x = k * radius
    th = np.asarray(pattern["theta_deg"], dtype=float)
    ib = int(np.argmin(np.abs(th - 180.0)))
    e0 = np.array([v if v is not None else 0.0 for v in pattern["series"][0]["E"]], dtype=float)
    denom = float(e0[ib] ** 2) or 1.0
    anchor = q_back * x * x / (4.0 * math.pi)
    return _bistatic_from_pattern(pattern, const=anchor / denom)


def _monostatic_sphere(k, radius, q_back):
    x = k * radius
    sig = q_back * x * x / (4.0 * math.pi)             # σ/λ² обратного рассеяния (изотропно)
    aspect = list(np.arange(0.0, 360.0 + 1e-9, 5.0))
    return {"aspect_deg": aspect, "sigma": [float(sig)] * len(aspect),
            "unit": "σ/λ²", "available": True, "isotropic": True,
            "note": "Сфера: ЭПР не зависит от ракурса (изотропно)."}


def _monostatic_cluster(scatterers, rad, gt, nmax):
    """Обратное рассеяние vs ракурс падения: развёртка k̂ в плоскости xz, поляризация ⊥
    плоскости (ŷ). T-матрицы и матрица взаимодействия зависят только от геометрии и k —
    собираются ОДИН раз, затем все ракурсы решаются разом (одна факторизация Mfull)."""
    from green_tensor import gmm, vswf
    k = rad["k"]
    P = len(scatterers)
    modes = vswf.mode_list(nmax)
    blk = 2 * len(modes)
    Ts = [gmm._full_t(s, k, nmax) for s in scatterers]                 # дорого, но 1 раз
    pos = [np.asarray(s.position, dtype=float) for s in scatterers]
    Mfull = np.eye(P * blk, dtype=complex)
    for p in range(P):
        for q in range(P):
            if p == q:
                continue
            # сборка идентична gmm.solve_cluster: вектор трансляции направлен от
            # нового центра (p) к старому (q), d = pos_q − pos_p; согласованность
            # с библиотекой закреплена тестом (studio/tests_studio.py)
            A, B, _ = vswf.translation_block_closed(pos[q] - pos[p], k, nmax, "out")
            G = np.block([[A, B], [B, A]])
            Mfull[p * blk:(p + 1) * blk, q * blk:(q + 1) * blk] = -G @ Ts[q]

    aspect = np.arange(0.0, 180.0 + 1e-9, MONO_STEP_DEG)
    khats = [(math.sin(math.radians(a)), 0.0, math.cos(math.radians(a))) for a in aspect]
    pol = (0.0, 1.0, 0.0)
    # все правые части (падения) сразу → одна LU-факторизация Mfull
    D = np.stack([np.concatenate([gmm.plane_wave_coeffs(k, kh, pol, pos[p], nmax)
                                  for p in range(P)]) for kh in khats], axis=1)
    Asol = np.linalg.solve(Mfull, D)                                  # (P·blk, N_aspect)
    R = FAR_R_FACTOR / k
    sig = []
    for j, kh in enumerate(khats):
        a = Asol[:, j].reshape(P, blk)
        c = np.array([Ts[p] @ a[p] for p in range(P)])
        pback = (R * np.array([-kh[0], -kh[1], -kh[2]])).reshape(1, 3)
        E = gmm.scattered_field(scatterers, c, k, pback)[0]
        f2 = float(np.sum(np.abs(E) ** 2)) * R * R                    # |f|²
        s = f2 * k * k / math.pi                                      # σ/λ²
        sig.append(float(s) if np.isfinite(s) else None)
    return {"aspect_deg": list(aspect), "sigma": _finite(sig), "unit": "σ/λ²", "available": True,
            "note": "Развёртка ракурса падения в плоскости xz; поляризация ⊥ плоскости."}


_RCS_2D_NOTE = ("Бесконечный цилиндр — 2D-задача: 3D-ЭПР неприменима "
                "(рассеяние характеризуется шириной σ₂D на единицу длины). "
                "См. вкладки «ДН» (азимутальная) и «Сечения» (Q на диаметр).")


# --------------------------------------------------------------------------- #
# Сечения и свип по k (только сфера — каноном 01_sphere)
# --------------------------------------------------------------------------- #
def _sphere_cross_sections(body, k, gt):
    radius = float(body.get("radius", 1.0))
    a_norm, eps, mu = _layers(body)
    cs = gt.SphereSolver(radius, eps, a_norm=a_norm, miy=mu).cross_sections(k)
    return {kk: (float(v) if np.isfinite(v) else None) for kk, v in cs.items()}


def _sphere_sweep(body, scene, gt):
    cfg = scene.get("sweep") or {}
    if not cfg.get("enable", True):
        return None
    radius = float(body.get("radius", 1.0))
    a_norm, eps, mu = _layers(body)
    k_cur = _radiation(scene)["k"]
    k_min = float(cfg.get("k0_min", 0.25))
    k_max = float(cfg.get("k0_max", max(2.0 * k_cur, 6.0)))
    npts = max(8, min(int(cfg.get("points", 160)), MAX_SWEEP_POINTS))
    if not (0 < k_min < k_max):
        raise ComputeError("Свип: требуется 0 < k0_min < k0_max.")
    grid = np.linspace(k_min, k_max, npts)
    out = {"k0": _finite(grid), "q_sca": [], "q_ext": [], "q_abs": [], "q_back": []}
    solver = gt.SphereSolver(radius, eps, a_norm=a_norm, miy=mu)
    for k in grid:
        try:
            cs = solver.cross_sections(float(k))
        except Exception:
            cs = {}
        for key in ("q_sca", "q_ext", "q_abs", "q_back"):
            v = cs.get(key)
            out[key].append(float(v) if v is not None and np.isfinite(v) else None)
    return out


# --------------------------------------------------------------------------- #
# Бесконечный радиально-слоистый цилиндр (строгий ТФГ, Дайлис–Шабунин 2017)
# --------------------------------------------------------------------------- #
def _cyl_inf_nmax(x):
    return int(math.ceil(abs(x) + 4.0 * abs(x) ** (1.0 / 3.0) + 2.0)) + 2


def _cyl_inf_pattern(solver, k, theta_deg, nmax, normal):
    """Азимутальная ДН |T(φ)| бесконечного цилиндра. Нормальное падение — раздельно
    TM/TE: T(φ)=a₀+2·Σₙ aₙ·cos(nφ). Косое — со- и кросс-поляризация (Σ по n=−N..N)."""
    phi = np.linspace(0.0, 2.0 * math.pi, N_THETA, endpoint=False)
    raw = []
    if normal:
        for mode, name in (("TM", "TM (E∥оси)"), ("TE", "TE (H∥оси)")):
            a = [complex(solver.coeff(k, n, mode=mode)) for n in range(nmax + 1)]
            T = a[0] + 2.0 * sum(a[n] * np.cos(n * phi) for n in range(1, nmax + 1))
            raw.append((name, np.abs(T)))
    else:
        theta = math.radians(theta_deg)
        Tco = np.zeros_like(phi, dtype=complex)
        Tcr = np.zeros_like(phi, dtype=complex)
        for n in range(-nmax, nmax + 1):
            ae, ah = solver.coeff_oblique(k, theta, n)
            Tco = Tco + complex(ae) * np.exp(1j * n * phi)
            Tcr = Tcr + complex(ah) * np.exp(1j * n * phi)
        raw.append(("Со-поляризация", np.abs(Tco)))
        raw.append(("Кросс-поляризация", np.abs(Tcr)))
    ref = max((float(np.max(a)) for _, a in raw if a.size), default=1.0) or 1.0
    series = [{"name": nm, "E": _finite(a), "dB": _finite(_db(a, ref))} for nm, a in raw]
    return {"theta_deg": _finite(np.degrees(phi)), "series": series,
            "polar": True, "plane": "φ"}


def _cyl_inf_sweep(body, scene, gt, mode):
    cfg = scene.get("sweep") or {}
    if not cfg.get("enable", True):
        return None
    R = float(body.get("radius", 0.5))
    a_norm, eps, mu = _layers(body)
    k_cur = _radiation(scene)["k"]
    k_min = float(cfg.get("k0_min", 0.25))
    k_max = float(cfg.get("k0_max", max(2.0 * k_cur, 6.0)))
    npts = max(8, min(int(cfg.get("points", 160)), MAX_SWEEP_POINTS))
    if not (0 < k_min < k_max):
        raise ComputeError("Свип: требуется 0 < k0_min < k0_max.")
    grid = np.linspace(k_min, k_max, npts)
    solver = gt.LayeredCylinderSolver(R, eps, mu=mu, a_norm=a_norm)
    out = {"k0": _finite(grid), "q_sca": [], "q_ext": [], "q_abs": []}
    for kk in grid:
        try:
            cs = solver.cross_sections(float(kk), mode=mode)
        except Exception:
            cs = {}
        for key in ("q_sca", "q_ext", "q_abs"):
            v = cs.get(key)
            out[key].append(float(v) if v is not None and np.isfinite(v) else None)
    return out


def _cylinder_inf(body, scene, rad, gt, *, do_pattern, do_sweep):
    """Строгий бесконечный радиально-слоистый цилиндр (ТФГ/эквивалентные линии,
    Дайлис–Шабунин 2017) через LayeredCylinderSolver. Эффективности Q на диаметр,
    азимутальная ДН, опц. спектр Q(k). Ближнего поля нет → heatmap=None."""
    R = float(body.get("radius", 0.5))
    a_norm, eps, mu = _layers(body)
    k = rad["k"]
    nmax = _cyl_inf_nmax(k * R)
    theta_deg = float(body.get("theta_deg", 90.0))
    mode = body.get("mode", "TM")
    if mode not in ("TM", "TE"):
        mode = "TM"
    normal = abs(theta_deg - 90.0) < 1e-6
    solver = gt.LayeredCylinderSolver(R, eps, mu=mu, a_norm=a_norm)

    if normal:
        cs = {}
        for m in ("TM", "TE"):
            d = solver.cross_sections(k, mode=m)
            for key in ("q_sca", "q_ext", "q_abs"):
                v = d.get(key)
                cs[f"{key}_{m}"] = float(v) if v is not None and np.isfinite(v) else None
    else:
        d = solver.cross_sections(k, theta=math.radians(theta_deg))
        cs = {key: (float(d[key]) if np.isfinite(d[key]) else None)
              for key in ("q_sca", "q_ext", "q_abs")}

    out = {"cross_sections": cs, "heatmap": None,
           "pattern": _cyl_inf_pattern(solver, k, theta_deg, nmax, normal) if do_pattern else None,
           "sweep": _cyl_inf_sweep(body, scene, gt, mode) if do_sweep else None}
    return out, nmax


# --------------------------------------------------------------------------- #
# Ближнепольная тепловая карта (любая геометрия) через Cluster.scattered_field
# --------------------------------------------------------------------------- #
def _heatmap_n(want):
    """Размер сетки near-field карты (точек на сторону) из запроса, с защитой пределами."""
    try:
        n = int(want.get("heatmap_res") or MAX_HEATMAP)
    except (TypeError, ValueError):
        n = MAX_HEATMAP
    return max(HEATMAP_RES_MIN, min(HEATMAP_RES_MAX, n))


def _heatmap_sphere(body, rad, extent, n=MAX_HEATMAP):
    """Полное ближнее поле сферы (однородной И слоистой) методом ТФГ/Ми
    (см. nearfield_sphere): падающее+рассеянное снаружи, внутреннее поле по слоям
    (трансфер-матрица — Z/Y пересчитываются по слоям), сверено с каноном 01_sphere
    по q_sca. Металлические слои (включая закрытое ядро) МАСКИРУЮТСЯ (поле не проникает).
    Возвращает |E| и доминирующую комплексную компоненту (re/im) для анимации напряжённости."""
    from .nearfield_sphere import field_sphere_layered, _wiscombe
    k = rad["k"]
    radius = float(body.get("radius", 1.0))
    a_norm, eps, mu = _layers(body)
    # Ряд ВСВФ обрезан на nmax≈Wiscombe(k·R) (по размеру сферы); вне r≈nmax/k он недостоверен
    # (на дальних точках нужен nmax≈k·r). При высоком k0 (напр. Люнеберг k0=6π) экстент 2.5·R
    # выходит за зону сходимости → артефакты в углах вместо |E|→1. Обрезаем кадр до зоны
    # сходимости (показываем сферу + ближнюю зону/фокус), маска r>r_valid — подстраховка.
    nmax_fld = max(_wiscombe(k * radius), 4)
    r_valid = nmax_fld / k
    extent = float(min(extent, 1.05 * r_valid))
    xs = np.linspace(-extent, extent, n)
    zs = np.linspace(-extent, extent, n)
    X, Z = np.meshgrid(xs, zs)
    pts = np.stack([X.ravel(), np.zeros(X.size), Z.ravel()], axis=-1)
    E, mask = field_sphere_layered(a_norm, eps, mu, k, radius, pts,
                                   khat=tuple(rad["khat"]), pol=tuple(rad["pol"]))
    absE = np.linalg.norm(E, axis=1)
    dom = int(np.argmax(np.abs(rad["pol"])))
    fld = E[:, dom]
    invalid = np.linalg.norm(pts, axis=1) > r_valid
    mask = mask | invalid
    absE[mask] = np.nan; fld[mask] = np.nan          # металл + вне зоны сходимости — не показываем
    absE = absE.reshape(n, n); fld = fld.reshape(n, n)
    return {
        "x": _finite(xs), "y": _finite(zs),
        "absE": [_finite(row) for row in absE],
        "reField": [_finite(row) for row in fld.real],
        "imField": [_finite(row) for row in fld.imag],
        "plane": "xz", "extent": extent, "internal": True, "res": n,
        "metal_masked": bool(mask.any()), "quantity": "total",
        "core": "01_sphere (ТФГ, внутр. поле по слоям)",
    }


def _heatmap(scatterers, rad, gt, nmax, extent, masks, n=MAX_HEATMAP, cf=None):
    """РАССЕЯННОЕ поле |E_расс| в плоскости xz (y=0) для кластеров/несфер. Точки ВНУТРИ
    описывающих сфер маскируются (NaN): разложение по исходящим h_n^(1) там не сходится.

    Показывается именно РАССЕЯННОЕ поле, а не полное: внутреннее поле для кластеров не
    считается, а полное |E_пад+E_расс| у слабых рассеивателей ≈ |E_пад|=1 (однородная
    плоская волна) и структура рассеяния не видна. Рассеянное поле выявляет лепестки/зону.

    cf — предрешённые рассеянные коэффициенты (если переданы — кластер не пересолвится).
    Для анимации напряжённости возвращается доминирующая комплексная компонента
    рассеянного поля — re/imField; |E_расс| — полная норма вектора.
    """
    k = rad["k"]
    xs = np.linspace(-extent, extent, n)
    zs = np.linspace(-extent, extent, n)
    X, Z = np.meshgrid(xs, zs)
    pts = np.stack([X.ravel(), np.zeros(X.size), Z.ravel()], axis=-1)   # плоскость y=0 (xz)
    inside = np.zeros(X.size, dtype=bool)
    # маскируем не только тело, но и БЛИЖНЮЮ зону у описывающей сферы: там исходящее
    # мультипольное разложение велико и на границе применимости (точечные «всплески» у
    # рёбер давили нормировку, карта казалась однородной). 1.25·R — надёжное пром. поле.
    for cx, cz, r in masks:
        inside |= (X.ravel() - cx) ** 2 + (Z.ravel() - cz) ** 2 < (r * 1.25) ** 2
    cl = gt.Cluster(scatterers)
    if cf is None:
        _, cf, _ = cl.solve(k, tuple(rad["khat"]), tuple(rad["pol"]), nmax)
    Esc = cl.scattered_field(cf, k, pts)                                # рассеянное (валидно вне тел)
    absE = np.linalg.norm(Esc, axis=1)
    dom = int(np.argmax(np.abs(rad["pol"])))                            # доминирующая компонента
    fld = Esc[:, dom]
    absE[inside] = np.nan; fld[inside] = np.nan                         # маска внутренней области
    absE = absE.reshape(n, n); fld = fld.reshape(n, n)
    return {
        "x": _finite(xs), "y": _finite(zs),
        "absE": [_finite(row) for row in absE],
        "reField": [_finite(row) for row in fld.real],
        "imField": [_finite(row) for row in fld.imag],
        "plane": "xz", "extent": extent, "internal": False, "res": n,
        "quantity": "scattered",
    }


def _mask_disk(body):
    """Пересечение описывающей сферы тела с плоскостью y=0 -> (cx, cz, r) или None."""
    pos = body.get("position", [0.0, 0.0, 0.0])
    cy = float(pos[1])
    R = float(body.get("radius", 1.0)) or 1.0
    if abs(cy) >= R:
        return None
    return (float(pos[0]), float(pos[2]), math.sqrt(R * R - cy * cy))


def _body_half_extent(b):
    """Полу-габарит тела (поперёк, вдоль z) — для авто-экстента карты (только сферы)."""
    r = float(b.get("radius", 1.0))
    return r, r            # сфера


def _heatmap_extent(bodies, margin=1.35):
    """Квадратный полу-экстент near-field карты по РЕАЛЬНЫМ габаритам тел + позиция.
    (Раньше бралось b['radius']=1.0 по умолчанию для сфероида ⇒ экстент раздувался и
    структура поля «терялась» в кадре с почти однородной падающей волной.)"""
    reach = 0.0
    for b in bodies:
        z = abs(float((b.get("position") or [0.0, 0.0, 0.0])[2]))
        lat, axz = _body_half_extent(b)
        reach = max(reach, z + axz, lat)
    return margin * (reach or 1.0)


# --------------------------------------------------------------------------- #
# Главная точка входа
# --------------------------------------------------------------------------- #
def compute(scene: dict) -> dict:
    gt = _load_library()
    bodies = scene.get("bodies") or []
    if not bodies:
        raise ComputeError("Сцена пуста: добавьте хотя бы одно тело.")
    rad = _radiation(scene)
    want = scene.get("compute") or {}
    do_pattern = want.get("pattern", True)
    do_heatmap = want.get("heatmap", True)
    do_sweep = want.get("sweep", False)   # свип медленный (канон 01_sphere ×N) — по запросу
    do_monostatic = want.get("monostatic", False)   # моностатич. ЭПР кластера — по запросу
    hm_n = _heatmap_n(want)               # разрешение тепловой карты (точек на сторону)
    warnings = []

    single_sphere = len(bodies) == 1 and bodies[0].get("type", "sphere") == "sphere"
    single_cyl_inf = len(bodies) == 1 and bodies[0].get("type") == "cylinder_inf"
    if len(bodies) > 1 and any(b.get("type") == "cylinder_inf" for b in bodies):
        raise ComputeError("Бесконечный цилиндр — 2D-задача; поддерживается только как одиночное тело.")

    out = {"ok": True, "meta": {"k": rad["k"], "n_bodies": len(bodies),
                                "single_sphere": single_sphere, "warnings": warnings}}

    if single_cyl_inf:
        body = bodies[0]
        res, nmax = _cylinder_inf(body, scene, rad, gt, do_pattern=do_pattern, do_sweep=do_sweep)
        out["cross_sections"] = res["cross_sections"]
        if do_pattern:
            out["pattern"] = res["pattern"]
        out["sweep"] = res["sweep"]
        out["heatmap"] = res["heatmap"]
        out["bistatic"] = {"available": False, "note": _RCS_2D_NOTE}
        out["monostatic"] = {"available": False, "note": _RCS_2D_NOTE}
        out["meta"]["core"] = "цилиндр (ТФГ, слоистый бесконечный)"
        out["meta"]["nmax"] = nmax
        out["meta"]["cylinder_infinite"] = True
        return out

    if single_sphere:
        body = bodies[0]
        radius = float(body.get("radius", 1.0))
        if _hemisphere_cfg(body):
            # Полусфера на PEC-экране: полупространственная задача, значима ДН.
            # Интегральные величины полной сферы (сечения, ЭПР, спектр, ближнее
            # поле плоской волны) к полупространству неприменимы — не выводятся.
            note = ("Полусфера на PEC-экране: сечения, ЭПР, спектр Q(k) и тепловая "
                    "карта для полупространственной задачи не выводятся — значима "
                    "диаграмма направленности.")
            if do_pattern:
                out["pattern"] = _hemisphere_pattern(body, rad, gt)
            out["cross_sections"] = None
            out["bistatic"] = {"available": False, "note": note}
            out["monostatic"] = {"available": False, "note": note}
            out["sweep"] = None
            out["heatmap"] = None
            out["meta"]["core"] = "01_sphere + зеркальный образ (полусфера на PEC-экране)"
            out["meta"]["warnings"].append(note)
            return out
        if rad["problem"] == "antenna":
            # Задача ИЗЛУЧЕНИЯ (источник на поверхности сферы): значима диаграмма
            # направленности. Сечения рассеяния, ЭПР, свип Q(k) и ближнее поле плоской
            # волны для этой задачи не определены.
            note = ("Задача «антенна» (источник на сфере): сечения рассеяния и ЭПР "
                    "не определены — выводится только диаграмма направленности.")
            if do_pattern:
                out["pattern"] = _sphere_pattern(body, rad, gt)
            out["cross_sections"] = None
            out["bistatic"] = {"available": False, "note": note}
            out["monostatic"] = {"available": False, "note": note}
            out["sweep"] = None
            out["heatmap"] = None
            out["meta"]["core"] = "01_sphere (canonical, antenna)"
            out["meta"]["warnings"].append(note)
            return out
        cs = _sphere_cross_sections(body, rad["k"], gt)
        out["cross_sections"] = cs
        q_back = cs.get("q_back") or 0.0
        if do_pattern:
            out["pattern"] = _sphere_pattern(body, rad, gt)
            if rad["polarization"] == "linear":
                out["bistatic"] = _bistatic_sphere(out["pattern"], rad["k"], radius, q_back)
            else:
                out["bistatic"] = {"available": False,
                                   "note": "Бистатическая ЭПР — для линейной поляризации."}
        out["monostatic"] = _monostatic_sphere(rad["k"], radius, q_back)
        out["sweep"] = _sphere_sweep(body, scene, gt) if do_sweep else None
        nmax = _nmax_for(bodies, rad["k"])
        if do_heatmap:
            # любая сфера (однородная/слоистая): полное поле вкл. внутри (ТФГ/Ми по слоям),
            # металл-слои маскируются
            out["heatmap"] = _heatmap_sphere(body, rad, 2.5 * radius, n=hm_n)
        out["meta"]["core"] = "01_sphere (canonical)"
        return out

    # Общий случай: кластер СФЕР через GMM (точная теорема сложения Крузана–Стейна).
    # Считаем ТОЛЬКО запрошенное (важно для opt-in моностатики).
    scatterers = [_build_scatterer(b, gt) for b in bodies]
    nmax = _nmax_for(bodies, rad["k"])
    out["meta"]["core"] = "GMM (кластер сфер)"
    out["meta"]["nmax"] = nmax
    do_cross = want.get("cross", True)

    # кластер решается ОДИН раз; рассеянные коэф. cf переиспользуются для сечений/ДН/карты
    cf = None
    try:
        if do_cross or do_pattern or do_heatmap:
            from green_tensor import gmm
            khat, pol = tuple(rad["khat"]), tuple(rad["pol"])
            cl = gt.Cluster(scatterers)
            _, cf, _ = cl.solve(rad["k"], khat, pol, nmax)
            if do_cross:
                cs = gmm.cross_sections(scatterers, cf, rad["k"], khat, pol)
                out["cross_sections"] = {kk: (float(v) if np.isfinite(v) else None)
                                         for kk, v in cs.items()}
            if do_pattern:
                out["pattern"] = _cluster_pattern(scatterers, rad, gt, nmax, cf=cf)
                out["bistatic"] = _bistatic_cluster(out["pattern"], rad["k"])
            if do_heatmap:
                ext = _heatmap_extent(bodies)
                masks = [m for m in (_mask_disk(b) for b in bodies) if m]
                out["heatmap"] = _heatmap(scatterers, rad, gt, nmax, extent=ext, masks=masks,
                                          n=hm_n, cf=cf)
        out["sweep"] = None
        out["monostatic"] = _monostatic_cluster(scatterers, rad, gt, nmax) if do_monostatic else None
    except ValueError as exc:
        # библиотека отклоняет физически недопустимые сцены (например, пересечение
        # описывающих сфер — теорема сложения неприменима); показываем причину в UI
        raise ComputeError(str(exc))
    return out
