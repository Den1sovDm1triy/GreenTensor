# SPDX-License-Identifier: MIT

"""Тесты модуля feed.py — перехват мощности облучателя раскрывом (eta_1).

Арбитры корректности:
  * канон диссертации: eta_1 элемента Гюйгенса на полусфере = 7/8 = 0.875
    (формула eq:eta1-spillover, расчёт coef_n.xmcd);
  * замкнутая форма для ЭГ: (8 - (1+cos theta0)^3)/8;
  * предел узкого конуса: eta_1 -> D*(1-cos theta0)/2 (интеграл измеряет КНД);
  * элементарный вибратор на полусфере: eta_1 = 1/2 (симметрия ДН).

Запуск:
    python3 tests/test_feed.py        # standalone, без pytest
    pytest tests/test_feed.py         # если установлен pytest
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from green_tensor.feed import (  # noqa: E402
    cone_half_angle,
    dipole_intensity,
    eta1_huygens_closed,
    feed_directivity,
    huygens_intensity,
    spillover_efficiency,
)


def test_huygens_hemisphere_canon():
    """ЭГ на поверхности линзы (theta0 = pi/2): eta_1 = 7/8 — канон диссертации."""
    eta = spillover_efficiency(huygens_intensity, np.pi / 2)
    assert abs(eta - 0.875) < 1e-10, eta


def test_huygens_closed_form():
    """Квадратура совпадает с замкнутой формой ЭГ на сетке углов, вкл. 5 град."""
    for deg in (1, 5, 10, 30, 60, 90, 120, 179):
        t0 = np.radians(deg)
        num = spillover_efficiency(huygens_intensity, t0)
        ref = eta1_huygens_closed(t0)
        assert abs(num - ref) < 1e-10, (deg, num, ref)


def test_dipole_hemisphere():
    """Элементарный вибратор: ДН симметрична относительно экватора, eta_1(pi/2)=1/2."""
    eta = spillover_efficiency(dipole_intensity, np.pi / 2)
    assert abs(eta - 0.5) < 1e-10, eta


def test_limits():
    """eta_1(0) = 0, eta_1(pi) = 1, монотонность по theta0."""
    assert spillover_efficiency(huygens_intensity, 0.0) == 0.0
    assert abs(spillover_efficiency(huygens_intensity, np.pi) - 1.0) < 1e-12
    grid = [spillover_efficiency(huygens_intensity, t)
            for t in np.linspace(0.0, np.pi, 37)]
    assert all(b >= a for a, b in zip(grid, grid[1:]))


def test_small_cone_measures_directivity():
    """Узкий конус: eta_1(theta0)/((1-cos theta0)/2) -> D (для ЭГ D = 3)."""
    D = feed_directivity(huygens_intensity)
    assert abs(D - 3.0) < 1e-6, D
    t0 = np.radians(1.0)
    ratio = spillover_efficiency(huygens_intensity, t0) / ((1 - np.cos(t0)) / 2)
    assert abs(ratio - D) < 1e-3, (ratio, D)


def test_tabulated_pattern_matches_callable():
    """Табличная ДН (theta_grid, |F|^2) даёт тот же eta_1, что и callable."""
    theta_grid = np.linspace(0.0, np.pi, 4001)
    f2_values = np.array([huygens_intensity(t) for t in theta_grid])
    for deg in (5, 45, 90):
        t0 = np.radians(deg)
        num = spillover_efficiency((theta_grid, f2_values), t0)
        ref = spillover_efficiency(huygens_intensity, t0)
        assert abs(num - ref) < 1e-6, (deg, num, ref)


def test_cone_geometry():
    """theta0 = arcsin(a/d): d = a -> pi/2 (полусфера); 5 град -> d/a ~ 11.47."""
    assert abs(cone_half_angle(1.0, 1.0) - np.pi / 2) < 1e-12
    d_over_a = 1.0 / np.sin(np.radians(5.0))
    assert abs(cone_half_angle(1.0, d_over_a) - np.radians(5.0)) < 1e-12


ALL_TESTS = (
    test_huygens_hemisphere_canon,
    test_huygens_closed_form,
    test_dipole_hemisphere,
    test_limits,
    test_small_cone_measures_directivity,
    test_tabulated_pattern_matches_callable,
    test_cone_geometry,
)


if __name__ == "__main__":
    failed = 0
    for t in ALL_TESTS:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as exc:
            failed += 1
            print(f"FAIL {t.__name__}: {exc}")
    sys.exit(1 if failed else 0)
