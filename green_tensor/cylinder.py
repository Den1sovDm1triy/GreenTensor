# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""cylinder — рассеяние на круговом цилиндре.

БЕСКОНЕЧНЫЙ круговой цилиндр (нормальное падение) решается ТОЧНО разделением
переменных в цилиндрических координатах (2D-ряд Ми, Bohren & Huffman гл. 8):
TM (E вдоль оси) и TE (H вдоль оси) расцеплены, коэффициенты через цилиндрические
Бессель/Ханкель целого порядка. Это каноническая аналитика (одна координатная
поверхность ρ=a).

КОНЕЧНЫЙ цилиндр — несепарабелен (две координатные поверхности + ребро), требует
сшивания мод / EBCM (см. GreenTensor_Theory.tex, раздел «Конечный цилиндр»). Эта
ветка не выдаётся непроверенной заглушкой: фините-решатель поднимает
NotImplementedError, а цилиндроподобные КОНЕЧНЫЕ тела для GMM собираются из
аналитического семейства через decompose.py (стопка сфер) либо аппроксимируются
вытянутым сфероидом (spheroid.py).

Конвенция e^{-iωt}; исходящая волна ~ H_n^{(1)}. x=k·a — параметр размера.
"""
from __future__ import annotations

import numpy as np
import scipy.special as sp


def _coeffs(m: complex, x: float, nmax: int, mode: str):
    """Коэффициенты рассеяния бесконечного цилиндра a_n, n=0..nmax (TM или TE)."""
    n = np.arange(0, nmax + 1)
    Jx, Jpx = sp.jv(n, x), sp.jvp(n, x)
    Jmx, Jpmx = sp.jv(n, m * x), sp.jvp(n, m * x)
    Hx, Hpx = sp.hankel1(n, x), sp.h1vp(n, x)
    if mode == "TM":      # E ∥ оси (case I, Bohren & Huffman)
        num = Jmx * Jpx - m * Jx * Jpmx
        den = Jmx * Hpx - m * Hx * Jpmx
    elif mode == "TE":    # H ∥ оси (case II, Bohren & Huffman)
        num = m * Jmx * Jpx - Jx * Jpmx
        den = m * Jmx * Hpx - Hx * Jpmx
    else:
        raise ValueError("mode must be 'TM' or 'TE'")
    return num / den


def cross_sections_infinite(m: complex, x: float, mode: str = "TM", nmax: int | None = None) -> dict:
    """Эффективности бесконечного цилиндра (норм. на диаметр 2a): Q_sca, Q_ext, Q_abs."""
    if nmax is None:
        nmax = int(np.ceil(abs(x) + 4.0 * abs(x) ** (1.0 / 3.0) + 2.0)) + 2
    a = _coeffs(m, x, nmax, mode)
    q_sca = (2.0 / x) * (abs(a[0]) ** 2 + 2.0 * np.sum(np.abs(a[1:]) ** 2))
    q_ext = (2.0 / x) * np.real(a[0] + 2.0 * np.sum(a[1:]))
    return {"q_sca": float(q_sca), "q_ext": float(q_ext), "q_abs": float(q_ext - q_sca)}


def _nmax_default(x: float) -> int:
    return int(np.ceil(abs(x) + 4.0 * abs(x) ** (1.0 / 3.0) + 2.0)) + 2


def _layer_basis(n: int, z: complex):
    """(J_n, H_n^{(1)}, J_n', H_n^{(1)'}) при (комплексном) аргументе z.

    Базис {J_n, H_n^{(1)}} устойчив для комплексного z (поглощающие слои); все
    функции scipy (jv/hankel1/jvp/h1vp) поддерживают комплексный аргумент.
    """
    return sp.jv(n, z), sp.hankel1(n, z), sp.jvp(n, z), sp.h1vp(n, z)


def layered_coeff(eps, mu, a_norm, x: float, n: int, mode: str = "TM") -> complex:
    """Коэффициент рассеяния a_n радиально-СЛОИСТОГО цилиндра (нормальное падение).

    Метод эквивалентных линий передачи (ТФГ, Дайлис–Шабунин 2017) при h=0: E- и
    H-линии расцеплены, поэтому каждая поляризация — независимая скалярная
    рекурсия амплитуд (A_i, B_i) поля E_z (TM) или H_z (TE) в базисе {J_n, H_n^{(1)}}.

    eps, mu : проницаемости слоёв изнутри наружу (комплексные допустимы);
    a_norm  : нормированные радиусы границ слоёв (внешний = 1), по возрастанию;
    x       : параметр размера на внешнем радиусе, x = k0·R;
    mode    : 'TM' (E∥ось) или 'TE' (H∥ось).
    """
    eps = [complex(e) for e in eps]
    mu = [complex(m) for m in mu]
    nL = len(eps)
    if len(mu) != nL or len(a_norm) != nL:
        raise ValueError("eps, mu, a_norm must have equal length")
    mrel = [np.sqrt(eps[i] * mu[i]) for i in range(nL)]   # относит. показатель слоя
    # множитель в условии непрерывности тангенциального H (TM) или E (TE):
    fac = [(mrel[i] / mu[i] if mode == "TM" else mrel[i] / eps[i]) for i in range(nL)]

    # ядро регулярно на оси ⇒ только J_n (амплитуда B при H_n^{(1)} равна нулю)
    A, B = 1.0 + 0j, 0.0 + 0j
    for j in range(nL - 1):                # внутренние границы r_j
        u = x * a_norm[j]
        Ji, Hi, Jpi, Hpi = _layer_basis(n, mrel[j] * u)
        Jo, Ho, Jpo, Hpo = _layer_basis(n, mrel[j + 1] * u)
        val = A * Ji + B * Hi                       # непрерывность E_z (TM)/H_z (TE)
        der = fac[j] * (A * Jpi + B * Hpi)          # непрерывность H_φ (TM)/E_φ (TE)
        a21, a22 = fac[j + 1] * Jpo, fac[j + 1] * Hpo
        det = Jo * a22 - Ho * a21
        A = (a22 * val - Ho * der) / det
        B = (-a21 * val + Jo * der) / det

    # внешняя граница u = x (a_norm[-1] = 1): сшивание с вакуумом
    JN, HN, JpN, HpN = _layer_basis(n, mrel[-1] * x)
    G = fac[-1] * (A * JpN + B * HpN) / (A * JN + B * HN)   # поверхностная проводимость
    Jx, Jpx = sp.jv(n, x), sp.jvp(n, x)
    Hx, Hpx = sp.hankel1(n, x), sp.h1vp(n, x)
    return (Jpx - G * Jx) / (Hpx - G * Hx)


def _field_matrix(n: int, eps_i: complex, mu_i: complex, g_i: complex,
                  u: float, cos_t: float) -> np.ndarray:
    """4×4 матрица [a,b,c,d] → (E_z, H_z, E_φ, H_φ) на радиусе arg = u·g_i.

    (a,b) — амплитуды J_n/H_n^{(1)} продольного E_z; (c,d) — то же для H_z.
    Связь E↔H в строках E_φ, H_φ идёт через член p = i·n·cosθ/u (∝ n·k_z),
    обнуляющийся при нормальном падении (cosθ=0) или n=0. Нормир. ед. ε0=μ0=1, ω=k0.
    """
    z = u * g_i
    J, H = sp.jv(n, z), sp.hankel1(n, z)
    Jp, Hp = sp.jvp(n, z), sp.h1vp(n, z)
    p = -1j * n * cos_t / u          # член связи E↔H (∝ n·k_z); знак — конвенция e^{-i k_z z}
    s = 1.0 / g_i ** 2
    M = np.zeros((4, 4), dtype=complex)
    M[0] = [J, H, 0.0, 0.0]                                              # E_z
    M[1] = [0.0, 0.0, J, H]                                             # H_z
    M[2] = [1j * s * p * J, 1j * s * p * H,
            -1j * s * mu_i * g_i * Jp, -1j * s * mu_i * g_i * Hp]        # E_φ
    M[3] = [1j * s * eps_i * g_i * Jp, 1j * s * eps_i * g_i * Hp,
            1j * s * p * J, 1j * s * p * H]                              # H_φ
    return M


def layered_coeff_oblique(eps, mu, a_norm, x: float, theta: float, n: int):
    """Рассеянные коэффициенты (A_E, A_H) слоистого цилиндра при КОСОМ падении.

    Метод ТФГ при h≠0 (Дайлис–Шабунин): E- и H-поляризации связаны на границах.
    4×4 матрица передачи на слой (непрерывность E_z, H_z, E_φ, H_φ); ядро регулярно
    (амплитуды H_n^{(1)} равны нулю). Падение — TM (вертикальная поляризация, E∥ось),
    угол ``theta`` от оси цилиндра (theta=π/2 — нормальное падение). Возвращает
    (A_E, A_H): со-поляризацию (рассеянный E_z) и кросс-поляризацию (рассеянный H_z).
    """
    eps = [complex(e) for e in eps]
    mu = [complex(m) for m in mu]
    nL = len(eps)
    cos_t = float(np.cos(theta))
    g = [np.sqrt(eps[i] * mu[i] - cos_t ** 2) for i in range(nL)]   # радиальное число / k0

    # пропагируем базис амплитуд ядра [a1, c1] через внутренние границы
    Vmat = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 1.0], [0.0, 0.0]], dtype=complex)
    for j in range(nL - 1):
        u = x * a_norm[j]
        Min = _field_matrix(n, eps[j], mu[j], g[j], u, cos_t)
        Mout = _field_matrix(n, eps[j + 1], mu[j + 1], g[j + 1], u, cos_t)
        Vmat = np.linalg.solve(Mout, Min @ Vmat)

    # внешняя граница u=x: сшивание с вакуумом (host: eps=mu=1, g0=sinθ)
    g0 = np.sqrt(1.0 - cos_t ** 2)
    MN = _field_matrix(n, eps[-1], mu[-1], g[-1], x, cos_t)
    Mh = _field_matrix(n, 1.0, 1.0, g0, x, cos_t)
    a_inc = float(np.sin(theta))                 # E_n^i = sinθ (как в арбитре Kavaklıoğlu)
    # неизвестные [a1, c1, A_E, A_H]: MN·(a1 V_a + c1 V_c) = Mh·[a_inc, A_E, 0, A_H]
    A = np.empty((4, 4), dtype=complex)
    A[:, 0] = MN @ Vmat[:, 0]
    A[:, 1] = MN @ Vmat[:, 1]
    A[:, 2] = -Mh[:, 1]                           # −(коэф. рассеянного E_z, H_n)
    A[:, 3] = -Mh[:, 3]                           # −(коэф. рассеянного H_z, H_n)
    rhs = a_inc * Mh[:, 0]
    sol = np.linalg.solve(A, rhs)
    return sol[2], sol[3]                          # (A_E со-пол, A_H кросс-пол)


def cross_sections_oblique(eps, mu, a_norm, x: float, theta: float,
                           nmax: int | None = None) -> dict:
    """Эффективности слоистого цилиндра при КОСОМ падении (TM/вертикальная поляризация).

    theta — угол падения от оси цилиндра (π/2 — нормальное). Учитывается и со-, и
    кросс-поляризация рассеянного поля. Нормировка как у Bohren & Huffman (на 2a, на
    падающую компоненту E_z = sinθ); экстинкция — по со-поляризации (опт. теорема).
    Проверено сохранением энергии (без потерь ⇒ Q_abs≈0) и переходом к
    :func:`cross_sections_layered` (mode='TM') при theta=π/2.
    """
    st = float(np.sin(theta))
    if not st > 1e-9:
        raise ValueError("theta слишком близко к оси (sinθ→0): задача вырождена")
    if nmax is None:
        nmax = _nmax_default(x)
    AE = np.empty(nmax + 1, dtype=complex)
    AH = np.empty(nmax + 1, dtype=complex)
    for n in range(nmax + 1):
        e, h = layered_coeff_oblique(eps, mu, a_norm, x, theta, n)
        AE[n], AH[n] = e / st, h / st         # коэффициенты на единичную падающую E_z
    q_ext = -(2.0 / x) * np.real(AE[0] + 2.0 * np.sum(AE[1:]))
    q_sca = (2.0 / x) * (abs(AE[0]) ** 2 + abs(AH[0]) ** 2
                         + 2.0 * np.sum(np.abs(AE[1:]) ** 2 + np.abs(AH[1:]) ** 2))
    return {"q_sca": float(q_sca), "q_ext": float(q_ext), "q_abs": float(q_ext - q_sca)}


def cross_sections_layered(eps, mu, a_norm, x: float, mode: str = "TM",
                           nmax: int | None = None) -> dict:
    """Эффективности слоистого цилиндра (норм. на диаметр): Q_sca, Q_ext, Q_abs.

    При одном слое и μ=1 совпадает с :func:`cross_sections_infinite`.
    """
    if nmax is None:
        nmax = _nmax_default(x)
    a = np.array([layered_coeff(eps, mu, a_norm, x, n, mode) for n in range(nmax + 1)],
                 dtype=complex)
    q_sca = (2.0 / x) * (abs(a[0]) ** 2 + 2.0 * np.sum(np.abs(a[1:]) ** 2))
    q_ext = (2.0 / x) * np.real(a[0] + 2.0 * np.sum(a[1:]))
    return {"q_sca": float(q_sca), "q_ext": float(q_ext), "q_abs": float(q_ext - q_sca)}


def finite(*args, **kwargs):
    """Конечный цилиндр (сшивание мод / EBCM) — не реализован.

    Несепарабелен (ρ=a и z=±L/2 + ребро); требует связанной модальной системы с
    краевым условием Мейкснера. Для GMM конечные цилиндры собираются из сфер
    (decompose.py) или приближаются вытянутым сфероидом (spheroid.py).
    """
    raise NotImplementedError(
        "Конечный цилиндр требует мод-матчинга/EBCM (см. GreenTensor_Theory.tex, "
        "раздел «Конечный цилиндр»). Используйте decompose.py или spheroid.py для GMM."
    )
