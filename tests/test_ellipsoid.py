"""Проверки квазистатического эллипсоида (green_tensor/ellipsoid.py).

  [1] Σ L_i = 1; сфера ⇒ L_i = 1/3;
  [2] сфероид: численные L против замкнутой формы (вытянутый);
  [3] предел сферы: рэлеевские Q_sca/Q_abs против независимой аналитики (analytic_mie);
  [4] дипольная T_N(1) сферы == рэлеевский предел коэффициента Ми a_1;
  [5] конфокальное покрытие с eps_core=eps_shell ⇒ однородный эллипсоид.

Запуск: python3 tests/test_ellipsoid.py | pytest tests/test_ellipsoid.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_ROOT, "green_tensor"))

from analytic_mie import q_sca, q_abs, mie_ab  # noqa: E402
import ellipsoid as el  # noqa: E402


def test_depolarization_sum_and_sphere():
    print("\n[1] ΣL=1 и сфера L=1/3:")
    L = el.depolarization_factors(1.0, 1.0, 1.0)
    assert abs(sum(L) - 1.0) < 1e-9
    assert all(abs(li - 1.0 / 3.0) < 1e-7 for li in L), f"сфера L={L}"
    La, Lb, Lc = el.depolarization_factors(3.0, 2.0, 1.0)
    print(f"    L(3,2,1)=({La:.4f},{Lb:.4f},{Lc:.4f}) sum={La+Lb+Lc:.6f}")
    assert abs(La + Lb + Lc - 1.0) < 1e-9 and La < Lb < Lc


def test_spheroid_closed_form():
    """Вытянутый сфероид (a>b=c): численный L_a против замкнутой формы."""
    print("\n[2] сфероид L_a: численно vs замкнуто:")
    a, b = 3.0, 1.0
    e = math.sqrt(1 - (b / a) ** 2)
    La_closed = (1 - e**2) / e**2 * (-1 + (1 / (2 * e)) * math.log((1 + e) / (1 - e)))
    La_num, Lb, Lc = el.depolarization_factors(a, b, b)
    print(f"    L_a числ={La_num:.6f} замкн={La_closed:.6f}  d={abs(La_num-La_closed):.1e}")
    assert abs(La_num - La_closed) < 1e-6
    assert abs(Lb - Lc) < 1e-9 and abs(2 * Lb - (1 - La_num)) < 1e-9


def test_sphere_rayleigh_vs_analytic():
    """Эллипсоид(a=b=c) в пределе Рэлея == аналитический Ми при малом x."""
    print("\n[3] сфера-эллипсоид: рэлеевские Q vs аналитика:")
    for eps in (2.25, 2.25 + 0.05j):
        for x in (0.02, 0.05):
            rad = 1.0
            k = x / rad
            alpha = el.polarizability_homogeneous(rad, rad, rad, eps)[0]
            cs = el.rayleigh_cross_sections(alpha, k)
            q_sca_el = cs["c_sca"] / (math.pi * rad**2)
            q_abs_el = cs["c_abs"] / (math.pi * rad**2)
            ref_sca = q_sca(np.sqrt(eps), x, 4)
            d_sca = abs(q_sca_el - ref_sca) / ref_sca
            print(f"    eps={eps!s:>12} x={x}: dQsca={d_sca:.1e}", end="")
            assert d_sca < 1e-2, f"Q_sca Рэлей разошёлся: {d_sca:.1e}"
            if abs(np.imag(eps)) > 0:
                ref_abs = q_abs(np.sqrt(eps), x, 4)
                d_abs = abs(q_abs_el - ref_abs) / abs(ref_abs)
                print(f" dQabs={d_abs:.1e}", end="")
                assert d_abs < 1e-2
            print()


def test_dipole_t_vs_mie():
    """Дипольная T_N(1) сферы из α == −a_1 (Ми) в пределе Рэлея."""
    print("\n[4] дипольная T_N(1) vs Ми a_1:")
    rad, eps = 1.0, 2.25
    for x in (0.02, 0.05):
        k = x / rad
        alpha = el.polarizability_homogeneous(rad, rad, rad, eps)[0]
        tN1 = el.dipole_t_scalar(alpha, k)
        _, a, _ = mie_ab(np.sqrt(eps), x, 4)
        d = abs(tN1 - (-a[0])) / abs(a[0])
        print(f"    x={x}: T_N(1)={tN1:.3e}  −a_1={-a[0]:.3e}  d={d:.1e}")
        assert d < 1e-2


def test_coated_reduces_to_homogeneous():
    print("\n[5] покрытие eps_core=eps_shell ⇒ однородный:")
    a, b, c = 2.0, 1.5, 1.0
    eps = 3.0 + 0.1j
    a_h = el.polarizability_homogeneous(a, b, c, eps)
    a_c = el.polarizability_coated_confocal(a, b, c, f=0.4, eps1=eps, eps2=eps)
    d = np.max(np.abs(a_h - a_c) / np.abs(a_h))
    print(f"    max rel diff = {d:.1e}")
    assert d < 1e-6


_TESTS = [test_depolarization_sum_and_sphere, test_spheroid_closed_form,
          test_sphere_rayleigh_vs_analytic, test_dipole_t_vs_mie,
          test_coated_reduces_to_homogeneous]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ ellipsoid проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
