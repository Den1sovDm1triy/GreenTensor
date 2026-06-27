# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Проверки интеграции КОНЕЧНОГО цилиндра в GMM через разложение в сферы.

green_tensor.decompose.decompose_cylinder и green_tensor.solvers.FiniteCylinderSolver:
конечный цилиндр аппроксимируется кластером непересекающихся сфер, который решается
сферическим движком GMM (gmm.solve_cluster). Прямой full-wave решатель отсутствует
(несепарабелен) и честно поднимает NotImplementedError.
"""
from __future__ import annotations

import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import green_tensor as gt  # noqa: E402
from green_tensor import decompose as dc  # noqa: E402


def test_indicator():
    print("\n[1] индикатор конечного цилиндра (внутри/снаружи):")
    ind = dc.cylinder_indicator([0, 0, 0], radius=1.0, half_length=2.0, axis=2)
    pts = np.array([[0, 0, 0], [0.9, 0, 1.9], [1.1, 0, 0], [0, 0, 2.1]], dtype=float)
    got = ind(pts).tolist()
    assert got == [True, True, False, False], got
    print(f"    {got} — ок")


def test_decompose_non_overlap_inside():
    print("\n[2] разложение: сферы не пересекаются и лежат внутри:")
    fc = gt.FiniteCylinderSolver([0, 0, 0], radius=2.0, half_length=4.0, eps=[2.25], axis=2)
    scat, centers, r = fc.decompose(spacing=1.0)
    assert len(scat) > 1
    assert dc.min_separation(centers) >= 2 * r - 1e-12       # строгое непересечение
    assert bool(np.all(fc.indicator()(centers)))             # центры внутри тела
    print(f"    сфер={len(scat)} radius={r:.3f} min_sep={dc.min_separation(centers):.3f} (2r={2*r:.3f}) — ок")


def test_decompose_feeds_gmm():
    print("\n[3] разложение → GMM (кластер решается, поле конечно):")
    fc = gt.FiniteCylinderSolver([0, 0, 0], radius=2.0, half_length=3.0, eps=[2.25])
    scat, _, _ = fc.decompose(spacing=1.0)
    a, c, d = gt.Cluster(scat).solve(k=0.7, khat=[1, 0, 0], pol=[0, 0, 1], nmax=3)
    assert np.all(np.isfinite(c)) and c.shape[0] == len(scat)
    cs = gt.Cluster(scat).cross_sections(k=0.7, khat=[1, 0, 0], pol=[0, 0, 1], nmax=3)
    assert all(np.isfinite(v) for v in cs.values())
    print(f"    сфер={len(scat)}, |c| конечно; C_ext={cs['c_ext']:.5f} — ок")


def test_solver_wrapper_and_notimpl():
    print("\n[4] FiniteCylinderSolver.decompose == decompose_cylinder; full_wave → NotImplementedError:")
    center, radius, hl, eps = [0.0, 0.0, 0.0], 2.0, 3.0, [2.25]
    fc = gt.FiniteCylinderSolver(center, radius, hl, eps, axis=2)
    s1, c1, r1 = fc.decompose(spacing=1.1)
    s2, c2, r2 = dc.decompose_cylinder(center, radius, hl, 1.1, eps, axis=2)
    assert len(s1) == len(s2) and np.allclose(c1, c2) and abs(r1 - r2) < 1e-15
    try:
        fc.full_wave()
        raised = False
    except NotImplementedError:
        raised = True
    assert raised, "full_wave должен поднимать NotImplementedError"
    print(f"    обёртка согласована ({len(s1)} сфер); full_wave — честный NotImplementedError — ок")


def _raises(fn):
    try:
        fn()
        return False
    except ValueError:
        return True


def test_metal_guard():
    print("\n[5] запрет разложения металла (паразитные полости); исключения:")
    C, R, HL = [0, 0, 0], 2.0, 4.0
    # (a) диэлектрик — разрешено
    assert len(gt.FiniteCylinderSolver(C, R, HL, [2.25]).decompose(1.0)[0]) > 1
    # (b) металл (Re eps<0) во внешнем слое — запрет
    assert _raises(lambda: gt.FiniteCylinderSolver(C, R, HL, [-10 + 1j]).decompose(1.0))
    # (c) металлическое ЯДРО под диэлектрической оболочкой — разрешено (внешний слой диэлектрик)
    assert len(gt.FiniteCylinderSolver(C, R, HL, [-10 + 1j, 2.25],
                                       a_norm=[0.6, 1.0]).decompose(1.0)[0]) > 1
    # (d) allow_metal=True — осознанный обход
    assert len(gt.FiniteCylinderSolver(C, R, HL, [-10 + 1j]).decompose(1.0, allow_metal=True)[0]) > 1
    # (e) лосси-проводник (Re>0, большой Im): без k проходит, с k (скин-слой<радиуса) — запрет
    fcs = gt.FiniteCylinderSolver(C, R, HL, [10 + 200j])
    assert len(fcs.decompose(1.0)[0]) > 1                       # Re>0, скин не проверен
    assert _raises(lambda: fcs.decompose(1.0, k=5.0))           # скин-слой ≪ радиуса
    # (f) тот же guard в конусе
    assert _raises(lambda: gt.ConeSolver([0, 0, 0], [0, 0, 1], 0.5, 8.0, [-10 + 1j]).decompose(2.0))
    # (g) низкоуровневый детектор
    assert dc.is_metal_layer(-3.0) and not dc.is_metal_layer(2.25)
    assert dc.is_metal_layer(10 + 200j, k=5.0, radius=0.45)
    print("    OK — диэлектрик/ядро-в-оболочке/allow_metal разрешены; металл-внешний и скин запрещены")


def test_fcc_lattice():
    print("\n[6] ГЦК-решётка: плотнее куба, непересекающаяся, выше достижимая ε (МГ):")
    inside = dc.sphere_indicator([0, 0, 0], 3.0)
    lo, hi = np.array([-3.0, -3, -3]), np.array([3.0, 3, 3])
    V = (4.0 / 3.0) * np.pi * 27.0
    Cc, rc = dc.pack_spheres(inside, lo, hi, 1.0, 0.45, lattice="cubic")
    Cf, rf = dc.pack_spheres(inside, lo, hi, 1.0, 0.45, lattice="fcc")
    fc = dc.coverage_fraction(Cc, rc, V)
    ff = dc.coverage_fraction(Cf, rf, V)
    assert dc.min_separation(Cf) >= 2 * rf - 1e-12        # строгое непересечение
    assert ff > fc                                         # ГЦК плотнее
    assert (1 + 2 * ff) / (1 - ff) > (1 + 2 * fc) / (1 - fc)   # выше eps_max для МГ
    # интеграция: FiniteCylinderSolver(lattice='fcc') даёт больше непересекающихся сфер
    sc, _, _ = gt.FiniteCylinderSolver([0, 0, 0], 2.0, 3.0, [2.0]).decompose(1.2, lattice="cubic")
    sf, cf2, rf2 = gt.FiniteCylinderSolver([0, 0, 0], 2.0, 3.0, [2.0]).decompose(1.2, lattice="fcc")
    assert len(sf) > len(sc) and dc.min_separation(cf2) >= 2 * rf2 - 1e-12
    print(f"    cubic f={fc:.3f} (ε≤{(1+2*fc)/(1-fc):.2f}) → fcc f={ff:.3f} (ε≤{(1+2*ff)/(1-ff):.2f}); "
          f"сфер {len(sc)}→{len(sf)} — ок")


if __name__ == "__main__":
    for fn in (test_indicator, test_decompose_non_overlap_inside,
               test_decompose_feeds_gmm, test_solver_wrapper_and_notimpl,
               test_metal_guard, test_fcc_lattice):
        fn()
    print("\nOK")
