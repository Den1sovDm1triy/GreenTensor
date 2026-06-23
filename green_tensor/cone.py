"""cone — рассеяние на конусе.

Полубесконечный круговой конус разделяется в сферо-конических координатах, но через
функции Лежандра НЕЦЕЛОЙ степени ν (корни из ГУ на поверхности θ=θ₀), а конечный
(усечённый) конус несепарабелен и требует сшивания мод (см. GreenTensor_Theory.tex,
раздел «Конечный конус»). Устойчивой библиотеки для P_ν^μ нецелой степени с нужной
точностью в scipy нет, поэтому прямой решатель НЕ выдаётся непроверенной заглушкой —
он поднимает NotImplementedError.

Практический и СТРОГО АНАЛИТИЧЕСКИЙ путь для конусоподобных тел в проекте: разложить
конус в кластер сфер (decompose.py) и решить через GMM — это и реализовано здесь
(cone_indicator + decompose_cone), и проверяется геометрически.
"""
from __future__ import annotations

import math

import numpy as np

from . import decompose as _dc


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
                   eps, fill: float = 0.45, a_norm=None, miy=None,
                   k: float | None = None, allow_metal: bool = False,
                   effective_medium: bool = False, lattice: str = "cubic"):
    """Разложить конус в непересекающиеся сферы и вернуть рассеиватели для GMM.

    Металл: разложение запрещено, если внешний слой металлический (см.
    :func:`decompose.reject_metal_packing`); k — для скин-слойной проверки,
    allow_metal=True — осознанный обход. effective_medium=True — поправка
    Максвелла–Гарнетта на ε (однородный диэлектрик; квазистатика)."""
    _dc.reject_metal_packing(eps, miy, _dc._lattice_radius(spacing, fill, lattice), k, allow_metal)
    apex = np.asarray(apex, dtype=float)
    ax = np.asarray(axis, dtype=float)
    ax = ax / np.linalg.norm(ax)
    R = height * math.tan(half_angle)
    # грубый bbox: вокруг отрезка apex..apex+height*ax, радиус R
    far = apex + height * ax
    lo = np.minimum(apex, far) - R
    hi = np.maximum(apex, far) + R
    inside = cone_indicator(apex, ax, half_angle, height)
    centers, radius = _dc.pack_spheres(inside, lo, hi, spacing, fill, lattice=lattice)
    eps_used = eps
    if effective_medium:
        body_volume = math.pi * R ** 2 * height / 3.0
        eps_used = _dc.effective_medium_eps(eps, miy, centers, radius, body_volume)
    return _dc.to_scatterers(centers, radius, eps_used, a_norm=a_norm, miy=miy), centers, radius


def full_wave(*args, **kwargs):
    """Прямой сферо-конический решатель (нецелые P_ν^μ) — не реализован.

    Требует устойчивого вычисления функций Лежандра нецелой степени и решения
    характеристического уравнения P_ν^μ(cosθ₀)=0 (semi-infinite) + сшивания мод
    (finite). См. GreenTensor_Theory.tex, раздел «Конечный конус». Для GMM используйте
    decompose_cone (строгая аналитика через сферы).
    """
    raise NotImplementedError(
        "Прямой решатель конуса требует функций Лежандра нецелой степени и мод-матчинга "
        "(см. GreenTensor_Theory.tex). Используйте decompose_cone для GMM."
    )
