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


if __name__ == "__main__":
    for fn in (test_indicator, test_decompose_non_overlap_inside,
               test_decompose_feeds_gmm, test_solver_wrapper_and_notimpl):
        fn()
    print("\nOK")
