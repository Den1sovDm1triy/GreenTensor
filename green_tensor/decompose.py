# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""decompose — разложение сложной геометрии в кластер аналитических примитивов.

Произвольное тело аппроксимируется набором НЕПЕРЕСЕКАЮЩИХСЯ сфер (описывающие сферы
не пересекаются — условие сходимости трансляционных рядов GMM, см.
GreenTensor_Theory.tex, раздел «Разложение»). Это даёт сборку сложной/«кирпичной»
формы из строго аналитического семейства без отдельного решателя для неё.

Сферы укладываются на кубическую решётку с шагом spacing и радиусом fill·spacing
(fill<0.5 ⇒ строгое непересечение). Возвращается список (центр, радиус), который
превращается в рассеиватели LayeredSphere для gmm.solve_cluster.
"""
from __future__ import annotations

import math

import numpy as np

from . import scatterer


def sphere_indicator(center, radius):
    c = np.asarray(center, dtype=float)
    return lambda P: np.sum((np.asarray(P, float) - c) ** 2, axis=-1) <= radius**2


def box_indicator(center, half):
    c = np.asarray(center, dtype=float)
    h = np.asarray(half, dtype=float)
    return lambda P: np.all(np.abs(np.asarray(P, float) - c) <= h, axis=-1)


def cylinder_indicator(center, radius, half_length, axis=2):
    c = np.asarray(center, dtype=float)
    perp = [i for i in range(3) if i != axis]
    def f(P):
        P = np.asarray(P, float) - c
        rad_ok = (P[..., perp[0]] ** 2 + P[..., perp[1]] ** 2) <= radius**2
        ax_ok = np.abs(P[..., axis]) <= half_length
        return rad_ok & ax_ok
    return f


def pack_spheres(inside_fn, bbox_min, bbox_max, spacing: float, fill: float = 0.45,
                 lattice: str = "cubic"):
    """Уложить непересекающиеся сферы внутрь тела (inside_fn).

    Возвращает (centers (M,3), radius). Сфера принимается, если её центр и 6 осевых
    точек на радиусе внутри тела. fill<0.5 ⇒ радиус меньше половины ближайшего соседа ⇒
    сферы строго не пересекаются (для любой решётки).

    lattice='cubic' — простая кубическая (шаг spacing, ближайший сосед = spacing,
    radius = fill·spacing). lattice='fcc' — гранецентрированная (spacing = сторона
    условной ячейки a, ближайший сосед = a/√2, radius = fill·a/√2): при тех же spacing
    и fill даёт примерно в √2 раз бо́льшую объёмную долю и лучшую конформность границе —
    важно для поправки эфф. среды (см. :func:`maxwell_garnett_eps`).
    """
    if not (0 < fill < 0.5):
        raise ValueError("fill must be in (0, 0.5) for non-overlapping spheres")
    bbox_min = np.asarray(bbox_min, float)
    bbox_max = np.asarray(bbox_max, float)

    if lattice == "cubic":
        nn = spacing
        axes = [np.arange(bbox_min[i] + spacing / 2, bbox_max[i], spacing) for i in range(3)]
        grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    elif lattice == "fcc":
        a = spacing
        nn = a / np.sqrt(2.0)                       # ближайший сосед в ГЦК
        cell_axes = [np.arange(bbox_min[i], bbox_max[i] + a, a) for i in range(3)]
        cells = np.stack(np.meshgrid(*cell_axes, indexing="ij"), axis=-1).reshape(-1, 3)
        basis = a * np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
        grid = (cells[:, None, :] + basis[None, :, :]).reshape(-1, 3)
        grid = grid[np.all((grid >= bbox_min) & (grid <= bbox_max), axis=1)]
    else:
        raise ValueError("lattice must be 'cubic' or 'fcc'")

    radius = fill * nn
    e = radius * np.eye(3)
    keep = inside_fn(grid)
    for i in range(3):
        keep &= inside_fn(grid + e[i]) & inside_fn(grid - e[i])
    return grid[keep], radius


def min_separation(centers) -> float:
    """Минимальное расстояние между центрами (для проверки непересечения)."""
    C = np.asarray(centers, float)
    if len(C) < 2:
        return np.inf
    best = np.inf
    for i in range(len(C)):
        dif = C[i + 1:] - C[i]
        if len(dif):
            best = min(best, np.min(np.sqrt(np.sum(dif**2, axis=1))))
    return best


def coverage_fraction(centers, radius: float, body_volume: float) -> float:
    """Доля объёма тела, покрытая сферами (Σ объёмов сфер / объём тела)."""
    vsph = (4.0 / 3.0) * np.pi * radius**3 * len(centers)
    return vsph / body_volume


def to_scatterers(centers, radius: float, eps, a_norm=None, miy=None):
    """Превратить упаковку в список рассеивателей LayeredSphere для GMM."""
    return [scatterer.LayeredSphere(c, radius, eps, a_norm=a_norm, miy=miy)
            for c in np.asarray(centers, float)]


# --------------------------------------------------------------------------- #
# Поправка эффективной среды (Максвелл–Гарнетт) + метрики упаковки
# Maxwell–Garnett effective-medium correction + packing metrics
# --------------------------------------------------------------------------- #
def maxwell_garnett_effective(eps_incl, f: float, eps_host: complex = 1.0) -> complex:
    """Эффективная ε смеси: включения eps_incl с объёмной долей f в среде eps_host.

    Формула Максвелла–Гарнетта (сферические включения):
        (eps_eff−eps_h)/(eps_eff+2eps_h) = f·(eps_i−eps_h)/(eps_i+2eps_h).
    """
    eps_i, eps_h = complex(eps_incl), complex(eps_host)
    beta = (eps_i - eps_h) / (eps_i + 2.0 * eps_h)
    return eps_h * (1.0 + 2.0 * f * beta) / (1.0 - f * beta)


def maxwell_garnett_eps(eps_eff, f: float, eps_host: complex = 1.0) -> complex:
    """Обратная задача МГ: ε включения, при котором эфф. среда (доля f) равна eps_eff.

    Используется для КОРРЕКЦИИ упаковки: сферы занимают долю f<1 объёма тела
    (остальное — вакуумные зазоры eps_host=1); подбираем ε сфер так, чтобы эффективная
    среда упаковки равнялась истинной ε тела. Тогда квазистатический отклик кластера
    совпадает со сплошным телом. Поднимает ValueError, если доли f не хватает для
    достижения eps_eff немагнитным диэлектриком (нужна более плотная упаковка).
    """
    if not 0.0 < f < 1.0:
        raise ValueError("объёмная доля f должна быть в (0, 1)")
    eps_e, eps_h = complex(eps_eff), complex(eps_host)
    L = (eps_e - eps_h) / (eps_e + 2.0 * eps_h)
    beta = L / f
    eps_i = eps_h * (1.0 + 2.0 * beta) / (1.0 - beta)
    if not (np.isfinite(eps_i) and eps_i.real > 0.0):
        raise ValueError(
            f"доля f={f:.3f} мала для достижения eps_eff={eps_eff} по Максвеллу–Гарнетту "
            "немагнитным диэлектриком: требуется более плотная упаковка (мельче spacing / "
            "плотнее решётка) либо результат non-physical (eps включения ≤ 0)."
        )
    return eps_i


def effective_medium_eps(eps, miy, centers, sphere_radius: float, body_volume: float):
    """ε сфер с поправкой Максвелла–Гарнетта под истинную ε ОДНОРОДНОГО тела.

    Возвращает список из одного элемента (для LayeredSphere). Применимо только к
    однородному немагнитному диэлектрику (один слой, mu=1): доля f=объём сфер/объём тела,
    eps_corr подбирается так, что эфф. среда упаковки = eps тела (см. maxwell_garnett_eps).
    """
    eps_list = list(eps)
    if len(eps_list) != 1:
        raise ValueError("поправка эфф. среды реализована только для однородного тела (один слой eps)")
    if miy is not None and any(abs(complex(m) - 1.0) > 1e-12 for m in miy):
        raise ValueError("поправка эфф. среды: магнитные слои (mu≠1) не поддержаны")
    f = coverage_fraction(centers, sphere_radius, body_volume)
    return [maxwell_garnett_eps(eps_list[0], f)]


def packing_report(centers, sphere_radius: float, body_volume: float,
                   nmax: int | None = None) -> dict:
    """Метрики упаковки: число сфер, объёмная доля, запас непересечения, оценка GMM.

    overlap_margin = min_sep/(2r) − 1 (>0 ⇒ строго не пересекаются; ~0 ⇒ касаются,
    сходимость GMM деградирует). При заданном nmax — размерность системы GMM и число
    элементов матрицы (оценка стоимости ~ (P·2K)²)."""
    centers = np.asarray(centers, dtype=float)
    n = len(centers)
    f = coverage_fraction(centers, sphere_radius, body_volume) if n else 0.0
    msep = min_separation(centers)
    margin = (msep / (2.0 * sphere_radius) - 1.0) if (n > 1 and sphere_radius > 0) else float("inf")
    rep = {"n_spheres": n, "filling": float(f), "min_sep": float(msep),
           "overlap_margin": float(margin), "sphere_radius": float(sphere_radius)}
    if nmax:
        dim = n * 2 * nmax * (nmax + 2)
        rep["gmm_dim"] = int(dim)
        rep["gmm_matrix_elems"] = int(dim) ** 2
    return rep


def is_metal_layer(eps, mu=1.0, *, k: float | None = None, radius: float | None = None) -> bool:
    """Является ли слой «металлическим» (поле в него не проникает).

    Критерий: Re(eps)<0 (плазмонный/металлический отклик) ИЛИ скин-слой меньше
    радиуса сферы δ=1/(k·Im√(εμ)) < radius (хороший проводник). Скин-слойная проверка
    требует k и radius; без них проверяется только Re(eps)<0.
    """
    eps = complex(eps)
    mu = complex(mu)
    if eps.real < 0.0:
        return True
    if k and radius and radius > 0:
        n_im = abs(np.sqrt(eps * mu).imag)
        if n_im > 0 and 1.0 / (k * n_im) < radius:
            return True
    return False


def reject_metal_packing(eps, miy, sphere_radius: float, k: float | None,
                         allow_metal: bool) -> None:
    """Запретить разложение в сферы, если ВНЕШНИЙ слой металлический (зазоры обнажены).

    Упаковка непересекающихся сфер оставляет вакуумные зазоры. Для металла это
    «губка»: зазоры — открытые полости, дающие паразитные переотражения и резонансы,
    которых у сплошного тела нет. Допустимо лишь металлическое ЯДРО, полностью
    закрытое диэлектрической внешней оболочкой (eps[-1] — диэлектрик): тогда зазоры
    касаются только диэлектрика, а металл — замкнутый объём внутри каждой сферы.
    """
    if allow_metal:
        return
    eps_list = list(eps)
    mu_outer = (list(miy)[-1] if miy is not None else 1.0)
    if is_metal_layer(eps_list[-1], mu_outer, k=k, radius=sphere_radius):
        raise ValueError(
            "Разложение в сферы недопустимо: внешний слой тела металлический "
            f"(eps_внеш={complex(eps_list[-1])}). Упаковка оставляет зазоры-полости, в "
            "которых возникают паразитные переотражения (поле не проникает в металл, "
            "но проникает в зазоры) — это металлическая «губка», а не сплошное тело. "
            "Допустимо: металлическое ЯДРО под диэлектрической внешней оболочкой "
            "(eps[-1] — диэлектрик), либо одиночное замкнутое тело "
            "(SphereSolver/LayeredSphere с металлическим слоем), либо поверхностный "
            "метод (EBCM/MoM). Осознанный обход (кластер металлических сфер): allow_metal=True."
        )


def _lattice_radius(spacing: float, fill: float, lattice: str) -> float:
    """Радиус сферы для решётки: cubic → fill·spacing, fcc → fill·spacing/√2."""
    return fill * spacing / (np.sqrt(2.0) if lattice == "fcc" else 1.0)


def decompose_cylinder(center, radius: float, half_length: float, spacing: float,
                       eps, *, axis: int = 2, fill: float = 0.45, a_norm=None, miy=None,
                       k: float | None = None, allow_metal: bool = False,
                       effective_medium: bool = False, lattice: str = "cubic"):
    """Разложить КОНЕЧНЫЙ круговой цилиндр в непересекающиеся сферы для GMM.

    center — центр; radius — радиус; half_length — половина длины вдоль оси; axis — ось
    цилиндра (0/1/2). Возвращает (scatterers, centers, sphere_radius). Сферический
    GMM (gmm.solve_cluster) собирает поле кластера — путь для конечного цилиндра в
    отсутствие прямого full-wave решателя (см. cylinder.finite).

    Металл: разложение запрещено, если внешний слой металлический (см.
    :func:`reject_metal_packing`); задайте k для скин-слойной проверки, allow_metal=True
    для осознанного обхода. effective_medium=True — поправка Максвелла–Гарнетта на ε
    (однородный диэлектрик; квазистатика), см. :func:`effective_medium_eps`."""
    reject_metal_packing(eps, miy, _lattice_radius(spacing, fill, lattice), k, allow_metal)
    center = np.asarray(center, dtype=float)
    half = np.full(3, radius, dtype=float)
    half[axis] = half_length
    lo, hi = center - half, center + half
    inside = cylinder_indicator(center, radius, half_length, axis)
    centers, sphere_radius = pack_spheres(inside, lo, hi, spacing, fill, lattice=lattice)
    eps_used = eps
    if effective_medium:
        body_volume = np.pi * radius ** 2 * (2.0 * half_length)
        eps_used = effective_medium_eps(eps, miy, centers, sphere_radius, body_volume)
    return to_scatterers(centers, sphere_radius, eps_used, a_norm=a_norm, miy=miy), centers, sphere_radius


def cone_indicator(apex, axis, half_angle: float, height: float):
    """Индикатор конечного кругового конуса: вершина apex, ось axis, полуугол, высота."""
    apex = np.asarray(apex, dtype=float)
    ax = np.asarray(axis, dtype=float)
    ax = ax / np.linalg.norm(ax)
    tan_a = math.tan(half_angle)

    def f(P):
        v = np.asarray(P, float) - apex
        axial = v @ ax
        perp = v - axial[..., None] * ax
        rad = np.sqrt(np.sum(perp**2, axis=-1))
        return (axial >= 0.0) & (axial <= height) & (rad <= axial * tan_a)

    return f


def decompose_cone(apex, axis, half_angle: float, height: float, spacing: float,
                   eps, *, fill: float = 0.45, a_norm=None, miy=None,
                   k: float | None = None, allow_metal: bool = False,
                   effective_medium: bool = False, lattice: str = "cubic"):
    """Разложить КОНЕЧНЫЙ круговой конус в непересекающиеся сферы для GMM.

    apex — вершина; axis — ось (от вершины к основанию); half_angle — полуугол; height —
    высота. Возвращает (scatterers, centers, sphere_radius). Это decompose-fallback
    конусоподобных тел; строгий полноволновой путь конуса — EBCM-примитив
    :class:`green_tensor.Cone` в :class:`green_tensor.solvers.Cluster`.

    Металл: разложение запрещено, если внешний слой металлический (см.
    :func:`reject_metal_packing`); k — для скин-слойной проверки, allow_metal=True —
    осознанный обход. effective_medium=True — поправка Максвелла–Гарнетта на ε
    (однородный диэлектрик; квазистатика), см. :func:`effective_medium_eps`."""
    reject_metal_packing(eps, miy, _lattice_radius(spacing, fill, lattice), k, allow_metal)
    apex = np.asarray(apex, dtype=float)
    ax = np.asarray(axis, dtype=float)
    ax = ax / np.linalg.norm(ax)
    R = height * math.tan(half_angle)
    far = apex + height * ax                      # грубый bbox: отрезок apex..far ± R
    lo = np.minimum(apex, far) - R
    hi = np.maximum(apex, far) + R
    inside = cone_indicator(apex, ax, half_angle, height)
    centers, radius = pack_spheres(inside, lo, hi, spacing, fill, lattice=lattice)
    eps_used = eps
    if effective_medium:
        body_volume = math.pi * R ** 2 * height / 3.0
        eps_used = effective_medium_eps(eps, miy, centers, radius, body_volume)
    return to_scatterers(centers, radius, eps_used, a_norm=a_norm, miy=miy), centers, radius
