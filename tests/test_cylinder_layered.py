# SPDX-License-Identifier: MIT

"""Проверки многослойного цилиндра (ТФГ, нормальное падение).

green_tensor/cylinder.py: layered_coeff / cross_sections_layered и обёртки
green_tensor.solvers.LayeredCylinderSolver / solve_layered_cylinder. Эталон —
независимый арбитр analytic_cylinder (Bohren & Huffman, однородный цилиндр).
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

import analytic_cylinder as arb  # noqa: E402
import green_tensor as gt  # noqa: E402
from green_tensor import cylinder as cyl  # noqa: E402


def _close(a, b, tol=1e-10):
    return abs(a - b) <= tol * (1.0 + abs(b))


def test_layered_n1_vs_arbiter():
    print("\n[1] слоистый N=1 == независимый арбитр B&H (TM/TE, поглощение, магн.):")
    for eps, mu in [(4.0, 1.0), (2.25 + 0.3j, 1.0), (4.0, 1.5)]:
        for x in (1.5, 3.0):
            for pol in ("TM", "TE"):
                la = cyl.cross_sections_layered([eps], [mu], [1.0], x, pol)
                ar = arb.cross_sections(eps, mu, x, pol)
                for key in ("q_sca", "q_ext", "q_abs"):
                    assert _close(la[key], ar[key]), (eps, mu, x, pol, key, la, ar)
    print("    OK — машинная точность")


def test_subdivision_invariance():
    print("\n[2] инвариантность к подразбиению (1 слой == 4 одинаковых):")
    eps, mu, x = 3.0 + 0.2j, 1.0, 2.5
    one = cyl.cross_sections_layered([eps], [mu], [1.0], x, "TM")
    four = cyl.cross_sections_layered([eps] * 4, [mu] * 4, [0.25, 0.5, 0.75, 1.0], x, "TM")
    for key in one:
        assert _close(one[key], four[key]), (key, one, four)
    print(f"    OK — Q_sca={one['q_sca']:.6f}")


def test_energy_conservation_lossless():
    print("\n[3] сохранение энергии, непоглощающий 3-слойный (Q_abs ~ 0):")
    for pol in ("TM", "TE"):
        cs = cyl.cross_sections_layered([2.0, 4.0, 1.5], [1.0, 1.0, 1.0],
                                        [0.4, 0.7, 1.0], 2.0, pol)
        assert abs(cs["q_abs"]) < 1e-10, (pol, cs)
    print("    OK — Q_abs ~ 1e-16")


def test_coated_absorption():
    print("\n[4] покрытый цилиндр (поглощающее ядро + прозрачная оболочка) Q_abs>0:")
    cs = cyl.cross_sections_layered([4.0 + 1.0j, 2.0], [1.0, 1.0], [0.6, 1.0], 2.0, "TM")
    assert cs["q_abs"] > 0 and cs["q_sca"] > 0, cs
    print(f"    OK — Q_sca={cs['q_sca']:.4f} Q_abs={cs['q_abs']:.4f}")


def test_solver_wrapper_matches():
    print("\n[5] LayeredCylinderSolver == cross_sections_layered == solve_layered_cylinder:")
    eps, mu, a_norm = [4.0 + 0.5j, 2.0], [1.0, 1.0], [0.6, 1.0]
    radius, k = 1.0, 2.0
    x = k * radius
    direct = cyl.cross_sections_layered(eps, mu, a_norm, x, "TM")
    oo = gt.LayeredCylinderSolver(radius, eps, mu=mu, a_norm=a_norm).cross_sections(k, mode="TM")
    fac = gt.solve_layered_cylinder(radius, eps, k, mu=mu, a_norm=a_norm, mode="TM")
    for key in ("q_sca", "q_ext", "q_abs"):
        assert _close(oo[key], direct[key]) and _close(fac[key], direct[key]), (key, oo, fac, direct)
    print(f"    OK — Q_ext={oo['q_ext']:.6f}")


def test_solver_single_layer_matches_homogeneous():
    print("\n[6] LayeredCylinderSolver(1 слой) == однородный CylinderSolver:")
    eps, radius, k = 4.0, 1.0, 2.0
    lay = gt.LayeredCylinderSolver(radius, [eps]).cross_sections(k, mode="TM")
    hom = gt.CylinderSolver(radius, eps).cross_sections(k, mode="TM")
    for key in ("q_sca", "q_ext", "q_abs"):
        assert _close(lay[key], hom[key]), (key, lay, hom)
    print(f"    OK — Q_sca={lay['q_sca']:.6f}")


def test_oblique_n1_vs_arbiter():
    print("\n[7] КОСОЕ падение, N=1 == арбитр Wait/Kavaklıoğlu (со- и кросс-пол):")
    worst = 0.0
    for eps in (4.0, 2.25, 3.0 + 0.4j):
        for x in (1.5, 3.0):
            for thd in (20.0, 45.0, 60.0, 80.0):
                th = np.deg2rad(thd)
                for n in range(0, 7):
                    AE, AH = cyl.layered_coeff_oblique([eps], [1.0], [1.0], x, th, n)
                    Aco, Acr = arb.oblique_coeffs_TM(eps, 1.0, x, th, n)
                    worst = max(worst, abs(AE - Aco), abs(AH - Acr))
    assert worst < 1e-9, worst
    print(f"    OK — max|Δ|={worst:.2e} (машинная точность)")


def test_oblique_reduces_to_normal():
    print("\n[8] косое → нормальное (θ=90°): кросс-пол → 0, со-пол == Фаза A:")
    wa = wr = 0.0
    for eps in (4.0, 3.0 + 0.4j):
        for x in (1.5, 3.0):
            for n in range(0, 6):
                AE, AH = cyl.layered_coeff_oblique([eps], [1.0], [1.0], x, np.pi / 2, n)
                an = cyl.layered_coeff([eps], [1.0], [1.0], x, n, "TM")
                wa = max(wa, abs(AH))
                wr = max(wr, abs(abs(AE) - abs(an)))
    assert wa < 1e-12 and wr < 1e-10, (wa, wr)
    print(f"    OK — max|A_H|={wa:.1e}, max||A_E|-|a_n||={wr:.1e}")


def test_oblique_subdivision_invariance():
    print("\n[9] косое падение: инвариантность к подразбиению (1 == 3 одинаковых):")
    eps, th, x = 3.0 + 0.2j, np.deg2rad(50.0), 2.5
    worst = 0.0
    for n in range(0, 5):
        one = cyl.layered_coeff_oblique([eps], [1.0], [1.0], x, th, n)
        three = cyl.layered_coeff_oblique([eps] * 3, [1.0] * 3, [1 / 3, 2 / 3, 1.0], x, th, n)
        worst = max(worst, abs(one[0] - three[0]), abs(one[1] - three[1]))
    assert worst < 1e-9, worst
    print(f"    OK — max|Δ|={worst:.2e}")


def test_oblique_energy_conservation():
    print("\n[10] косое падение: сохранение энергии (непоглощ. многослойный, Q_abs~0):")
    for thd in (30.0, 45.0, 70.0):
        cs = cyl.cross_sections_oblique([2.0, 4.0, 1.5], [1.0, 1.0, 1.0],
                                        [0.4, 0.7, 1.0], 2.0, np.deg2rad(thd))
        assert abs(cs["q_abs"]) < 1e-9, (thd, cs)
    print("    OK — Q_abs ~ 1e-15 на всех углах")


def test_oblique_cross_sections_reduce_and_wrapper():
    print("\n[11] косое: cross_sections(θ=90°) == Фаза A; обёртки солвера/фасада:")
    eps, mu, a_norm = [4.0 + 0.3j, 2.0], [1.0, 1.0], [0.6, 1.0]
    radius, k = 1.0, 2.0
    x = k * radius
    # при θ=90° косая формула обязана совпасть с нормальной TM
    obl90 = cyl.cross_sections_oblique(eps, mu, a_norm, x, np.pi / 2)
    norm = cyl.cross_sections_layered(eps, mu, a_norm, x, "TM")
    for key in ("q_sca", "q_ext", "q_abs"):
        assert _close(obl90[key], norm[key], tol=1e-9), (key, obl90, norm)
    # обёртки: LayeredCylinderSolver(theta) и solve_layered_cylinder(theta)
    th = np.deg2rad(40.0)
    direct = cyl.cross_sections_oblique(eps, mu, a_norm, x, th)
    oo = gt.LayeredCylinderSolver(radius, eps, mu=mu, a_norm=a_norm).cross_sections(k, theta=th)
    fac = gt.solve_layered_cylinder(radius, eps, k, mu=mu, a_norm=a_norm, theta=th)
    for key in ("q_sca", "q_ext", "q_abs"):
        assert _close(oo[key], direct[key]) and _close(fac[key], direct[key]), (key, oo, fac, direct)
    print(f"    OK — θ=90° совпало с Фазой A; обёртки согласованы (Q_ext(40°)={direct['q_ext']:.4f})")


if __name__ == "__main__":
    arb.__dict__  # noqa: B018  (импорт-проверка)
    for fn in (test_layered_n1_vs_arbiter, test_subdivision_invariance,
               test_energy_conservation_lossless, test_coated_absorption,
               test_solver_wrapper_matches, test_solver_single_layer_matches_homogeneous,
               test_oblique_n1_vs_arbiter, test_oblique_reduces_to_normal,
               test_oblique_subdivision_invariance, test_oblique_energy_conservation,
               test_oblique_cross_sections_reduce_and_wrapper):
        fn()
    print("\nOK")
