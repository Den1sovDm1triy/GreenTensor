# SPDX-License-Identifier: MIT

"""Абсолютная бистатическая ЭПР σ(θ) — проверка амплитудного префактора.

Ядро строит ненормированные амплитуды S_1, S_2 и абсолютную σ(θ) = (4π/k²)|S|².
Единственная внутренняя проверка (без внешних таблиц): в точке обратного рассеяния
абсолютная бистатическая ЭПР равна уже валидированной интегральной величине —

    σ_E(θ=π) = σ_H(θ=π) = π R² · Q_back,

где Q_back проверен против ряда Ми в других тестах (test_sphere_multilayer, test_tmatrix).
Тем самым фиксируется корректность множителя 4π/k² для однородной, PEC, магнитодиэлект-
рической и многослойной сфер.

Запуск:
    python3 tests/test_bistatic_rcs.py     # standalone
    pytest tests/test_bistatic_rcs.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from green_tensor import bistatic_rcs_sphere, far_field_sphere, solve_sphere  # noqa: E402

# (радиус R, eps слоёв, a_norm, miy, k) — изнутри наружу
_CASES = [
    ("однородная eps=4",        1.0,   [4.0],                 None,         None,   3.0),
    ("однородная eps=4 (x=8)",  1.0,   [4.0],                 None,         None,   8.0),
    ("PEC",                     1.0,   [1e7],                 None,         None,   5.0),
    ("магнитодиэлектрик",       1.0,   [4.0],                 None,         [2.0],  2.0),
    ("двухслойная [4,2]",       1.0,   [4.0, 2.0],            [0.5, 1.0],   None,   4.0),
    ("PEC + укрытие Rogers",    3.587e-3, [1e7, 10.2 + 0.02346j], [0.5 / 0.75, 1.0], None, 0.75 / 3.587e-3),
    ("линза Люнеберга k0a=15",  1.0, [1.96, 1.84, 1.64, 1.36, 1.0],
     [0.2, 0.4, 0.6, 0.8, 1.0], None, 15.0),
]


def _backscatter_matches_qback(R, eps, a_norm, miy, k, rtol=1e-9):
    q_back = solve_sphere(radius=R, eps=eps, k=k, a_norm=a_norm, miy=miy)["q_back"]
    s = bistatic_rcs_sphere(R, eps, k, np.array([math.pi]), a_norm=a_norm, miy=miy)
    target = math.pi * R**2 * q_back
    ok_e = math.isclose(float(s["sigma_E"][0]), target, rel_tol=rtol)
    ok_h = math.isclose(float(s["sigma_H"][0]), target, rel_tol=rtol)
    return ok_e and ok_h, float(s["sigma_E"][0]), float(s["sigma_H"][0]), target


def test_backscatter_equals_pi_r2_qback():
    """σ_E(π) = σ_H(π) = π R² Q_back для всех классов сфер."""
    for name, R, eps, a_norm, miy, k in _CASES:
        ok, se, sh, tgt = _backscatter_matches_qback(R, eps, a_norm, miy, k)
        assert ok, f"{name}: sigma_E={se:.6e} sigma_H={sh:.6e} != piR^2 Qback={tgt:.6e}"


def test_far_field_consistent_with_rcs():
    """σ = 4π|rE|² при E0=1 В/м (rE = |S|/k) — связь поля и ЭПР."""
    R, eps, k = 1.0, [4.0], 3.0
    theta = np.linspace(0.01, math.pi, 64)
    s = bistatic_rcs_sphere(R, eps, k, theta)
    f = far_field_sphere(R, eps, k, theta, e0=1.0)
    np.testing.assert_allclose(4 * math.pi * f["rE_E"] ** 2, s["sigma_E"], rtol=1e-9)
    np.testing.assert_allclose(4 * math.pi * f["rE_H"] ** 2, s["sigma_H"], rtol=1e-9)


def test_field_scales_with_distance():
    """E(r) = rE / r."""
    R, eps, k, r = 1.0, [2.25], 5.0, 10.0
    theta = np.linspace(0.01, math.pi, 32)
    f = far_field_sphere(R, eps, k, theta, e0=1.0, r=r)
    np.testing.assert_allclose(f["E_E"] * r, f["rE_E"], rtol=1e-12)


if __name__ == "__main__":
    print("Абсолютная бистатическая ЭПР: σ_E(π)=σ_H(π)=π R² Q_back\n")
    all_ok = True
    for name, R, eps, a_norm, miy, k in _CASES:
        ok, se, sh, tgt = _backscatter_matches_qback(R, eps, a_norm, miy, k)
        all_ok &= ok
        print(f"  [{'OK' if ok else 'FAIL'}] {name:26s} σ_E={se:.6e}  π R² Q_back={tgt:.6e}  "
              f"отн.разн={abs(se - tgt) / tgt:.1e}")
    print("\nСвязь поля и ЭПР (σ = 4π|rE|²): ", end="")
    try:
        test_far_field_consistent_with_rcs()
        test_field_scales_with_distance()
        print("OK")
    except AssertionError as exc:
        all_ok = False
        print(f"FAIL — {exc}")
    print("\nИТОГ:", "все проверки пройдены" if all_ok else "ЕСТЬ ОШИБКИ")
    sys.exit(0 if all_ok else 1)
