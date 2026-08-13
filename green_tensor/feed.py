# SPDX-License-Identifier: MIT
"""GreenTensor — расчёт перехвата мощности облучателя раскрывом (spill-over).

Коэффициент использования поля облучателя eta_1 — доля мощности первичного
источника, попадающая в конус с полууглом theta0, под которым раскрыв (линза)
виден из фазового центра облучателя:

    eta_1(theta0) = int_0^theta0 |F(theta)|^2 sin(theta) d(theta)
                    / int_0^pi |F(theta)|^2 sin(theta) d(theta)

(осесимметричная ДН; общий случай отличается лишь интегрированием по phi).
Частный случай theta0 = pi/2 — облучатель на поверхности линзы (полусфера) —
соответствует формуле (eq:eta1-spillover) главы о ТФГ диссертации и расчёту
coef_n.xmcd из серии «Эффективность облучателя» (_CalculateReferences/).

Полезные пределы, используемые как внутренние тесты корректности:
  * eta_1(0) = 0, eta_1(pi) = 1;
  * eta_1(theta0 -> 0) -> D * (1 - cos theta0)/2, где D — КНД облучателя,
    то есть интеграл по узкому конусу напрямую измеряет КНД;
  * для элемента Гюйгенса замкнутая форма eta_1 = (8 - (1 + cos theta0)^3)/8
    даёт канонические 7/8 = 0.875 при theta0 = pi/2.

Feed spill-over module: fraction of the primary-feed power intercepted by an
aperture (lens) subtending half-angle theta0 at the feed phase centre.
Special case theta0 = pi/2 (feed on the lens surface) reproduces the
dissertation's eta_1 values (Huygens element: 0.875).
"""
from __future__ import annotations

from typing import Callable, Sequence, Union

import numpy as np
from scipy.integrate import quad

__all__ = [
    "huygens_intensity",
    "dipole_intensity",
    "spillover_efficiency",
    "eta1_huygens_closed",
    "feed_directivity",
    "cone_half_angle",
]

# Осесимметричная нормированная ДН по мощности |F(theta)|^2
IntensityPattern = Union[Callable[[float], float], tuple]


def huygens_intensity(theta: float) -> float:
    """|F(theta)|^2 элемента Гюйгенса (апертурный излучатель): ((1+cos)/2)^2."""
    return float(((1.0 + np.cos(theta)) / 2.0) ** 2)


def dipole_intensity(theta: float) -> float:
    """|F(theta)|^2 элементарного электрического вибратора: sin^2(theta)."""
    return float(np.sin(theta) ** 2)


def _as_callable(pattern: IntensityPattern) -> Callable[[float], float]:
    """Приводит ДН к вызываемому виду.

    Принимает либо callable |F(theta)|^2, либо пару (theta_grid, f2_values)
    с табличной ДН (линейная интерполяция; сетка должна покрывать [0, pi]).
    """
    if callable(pattern):
        return pattern
    theta_grid, f2_values = pattern
    theta_grid = np.asarray(theta_grid, dtype=float)
    f2_values = np.asarray(f2_values, dtype=float)
    if theta_grid[0] > 0.0 or theta_grid[-1] < np.pi:
        raise ValueError("табличная ДН должна покрывать отрезок [0, pi]")
    return lambda t: float(np.interp(t, theta_grid, f2_values))


def spillover_efficiency(pattern: IntensityPattern, theta0: float) -> float:
    """eta_1(theta0) — доля мощности облучателя в конусе полуугла theta0.

    pattern — ДН облучателя в СВОБОДНОМ пространстве (не системы с линзой):
    callable |F(theta)|^2 или пара (theta_grid, f2_values).
    theta0 — полуугол конуса раскрыва, рад, 0 <= theta0 <= pi.
    """
    if not 0.0 <= theta0 <= np.pi:
        raise ValueError("theta0 должен лежать в [0, pi]")
    f2 = _as_callable(pattern)
    integrand = lambda t: f2(t) * np.sin(t)
    total, _ = quad(integrand, 0.0, np.pi, limit=200)
    if total <= 0.0:
        raise ValueError("нулевая полная мощность ДН")
    partial, _ = quad(integrand, 0.0, theta0, limit=200)
    return partial / total


def eta1_huygens_closed(theta0: float) -> float:
    """Замкнутая форма eta_1 для элемента Гюйгенса: (8 - (1+cos theta0)^3)/8."""
    return (8.0 - (1.0 + np.cos(theta0)) ** 3) / 8.0


def feed_directivity(pattern: IntensityPattern) -> float:
    """КНД осесимметричного облучателя: D = 2*max|F|^2 / int |F|^2 sin d(theta)."""
    f2 = _as_callable(pattern)
    total, _ = quad(lambda t: f2(t) * np.sin(t), 0.0, np.pi, limit=200)
    peak = max(f2(t) for t in np.linspace(0.0, np.pi, 1801))
    return 2.0 * peak / total


def cone_half_angle(lens_radius: float, feed_distance: float) -> float:
    """Полуугол theta0 = arcsin(a/d), под которым линза радиуса a видна
    из фазового центра облучателя на расстоянии d >= a от центра линзы.
    d = a (облучатель на поверхности) даёт pi/2 — полусферу диссертации."""
    if feed_distance < lens_radius:
        raise ValueError("облучатель внутри линзы: d < a")
    return float(np.arcsin(lens_radius / feed_distance))
