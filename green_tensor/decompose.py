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


def pack_spheres(inside_fn, bbox_min, bbox_max, spacing: float, fill: float = 0.45):
    """Уложить непересекающиеся сферы внутрь тела (inside_fn) на кубической решётке.

    Возвращает (centers (M,3), radius). radius=fill·spacing, fill<0.5 ⇒ строго не пересекаются.
    Сфера принимается, если её центр и 6 осевых точек на радиусе внутри тела.
    """
    if not (0 < fill < 0.5):
        raise ValueError("fill must be in (0, 0.5) for non-overlapping spheres")
    radius = fill * spacing
    bbox_min = np.asarray(bbox_min, float)
    bbox_max = np.asarray(bbox_max, float)
    axes = [np.arange(bbox_min[i] + spacing / 2, bbox_max[i], spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
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


def decompose_cylinder(center, radius: float, half_length: float, spacing: float,
                       eps, *, axis: int = 2, fill: float = 0.45, a_norm=None, miy=None):
    """Разложить КОНЕЧНЫЙ круговой цилиндр в непересекающиеся сферы для GMM.

    center — центр; radius — радиус; half_length — половина длины вдоль оси; axis — ось
    цилиндра (0/1/2). Возвращает (scatterers, centers, sphere_radius). Сферический
    GMM (gmm.solve_cluster) собирает поле кластера — путь для конечного цилиндра в
    отсутствие прямого full-wave решателя (см. cylinder.finite)."""
    center = np.asarray(center, dtype=float)
    half = np.full(3, radius, dtype=float)
    half[axis] = half_length
    lo, hi = center - half, center + half
    inside = cylinder_indicator(center, radius, half_length, axis)
    centers, sphere_radius = pack_spheres(inside, lo, hi, spacing, fill)
    return to_scatterers(centers, sphere_radius, eps, a_norm=a_norm, miy=miy), centers, sphere_radius
