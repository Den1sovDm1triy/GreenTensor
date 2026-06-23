"""Проверки сфероида (green_tensor/spheroid.py, квазистатика).

  [1] замкнутая форма L_ax (prolate и oblate) == численный интеграл;
  [2] предел сферы (c_ax=a_eq) ⇒ L=1/3;
  [3] почти-сфера: дипольная T_N(1) ≈ Ми a_1 при малом x;
  [4] полноволновая ветка честно поднимает NotImplementedError.

Запуск: python3 tests/test_spheroid.py | pytest tests/test_spheroid.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from analytic_mie import mie_ab  # noqa: E402
from green_tensor import spheroid as sph  # noqa: E402


def test_closed_form_vs_numerical():
    print("\n[1] замкнутая L_ax vs численная:")
    for (a_eq, c_ax, name) in [(1.0, 3.0, "prolate"), (2.0, 0.7, "oblate")]:
        L_eq_n, _, L_ax_n = sph.depolarization(a_eq, c_ax)
        L_eq_c, L_ax_c = sph.depolarization_closed(a_eq, c_ax)
        print(f"    {name}: L_ax числ={L_ax_n:.6f} замкн={L_ax_c:.6f} d={abs(L_ax_n-L_ax_c):.1e}")
        assert abs(L_ax_n - L_ax_c) < 1e-7 and abs(L_eq_n - L_eq_c) < 1e-7
        assert abs(2 * L_eq_n + L_ax_n - 1.0) < 1e-9


def test_sphere_limit():
    print("\n[2] предел сферы:")
    L_eq, L_ax = sph.depolarization_closed(1.0, 1.0)
    assert abs(L_eq - 1.0 / 3.0) < 1e-12 and abs(L_ax - 1.0 / 3.0) < 1e-12
    Le, _, La = sph.depolarization(1.0, 1.0)
    assert abs(Le - 1.0 / 3.0) < 1e-7 and abs(La - 1.0 / 3.0) < 1e-7
    print("    c_ax=a_eq ⇒ L=1/3 — ок")


def test_near_sphere_dipole_vs_mie():
    print("\n[3] почти-сфера: дипольная T_N(1) vs Ми a_1:")
    a = 1.0
    eps = 2.25
    for x in (0.02, 0.05):
        k = x / a
        alpha = sph.polarizability(a, a, eps)[0]   # сфера
        tN1 = sph.dipole_t_scalar(alpha, k)
        _, am, _ = mie_ab(math.sqrt(eps), x, 4)
        d = abs(tN1 - (-am[0])) / abs(am[0])
        print(f"    x={x}: d={d:.1e}")
        assert d < 1e-2


def test_full_wave_not_implemented():
    print("\n[4] полноволновая ветка → NotImplementedError:")
    try:
        sph.full_wave()
        raised = False
    except NotImplementedError:
        raised = True
    assert raised, "full_wave должна честно поднимать NotImplementedError"
    print("    честный NotImplementedError — ок")


_TESTS = [test_closed_form_vs_numerical, test_sphere_limit,
          test_near_sphere_dipole_vs_mie, test_full_wave_not_implemented]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ spheroid проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
