"""Проверки публичного API green_tensor/solvers.py.

Каждый класс-решатель и функция-обёртка ``solve_*`` НЕ содержат своей математики —
тесты подтверждают, что они в точности воспроизводят уже проверенные нижележащие
функции (никакого расхождения), и что полноволновые ветви честно поднимают
NotImplementedError.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

import green_tensor as gt  # noqa: E402
from green_tensor import cylinder, ellipsoid, gmm, mie_core, spheroid, tmatrix, vswf  # noqa: E402


def _close(a, b, tol=1e-12):
    return abs(a - b) <= tol * (1.0 + abs(b))


def test_sphere_solver_matches_core():
    print("\n[1] SphereSolver == mie_core == solve_sphere:")
    radius, eps, k = 1.0, [2.25], 3.0
    sp = gt.SphereSolver(radius, eps)
    cs = sp.cross_sections(k)
    ref = mie_core.MieSphere(k0=k * radius, a=[1.0], eps=eps).cross_sections()
    facade = gt.solve_sphere(radius, eps, k)
    for key in ("q_sca", "q_ext", "q_abs", "q_back"):
        assert _close(cs[key], ref[key]), key
        assert _close(facade[key], ref[key]), key
    print(f"    Q_sca={cs['q_sca']:.6f} (== core, == facade) — ок")


def test_sphere_tmatrix_consistency():
    print("\n[2] SphereSolver.t_matrix.q_sca == cross_sections:")
    sp = gt.SphereSolver(1.0, [2.25])
    cs = sp.cross_sections(3.0)
    tm = sp.t_matrix(3.0)
    assert isinstance(tm, tmatrix.DiagonalTMatrix)
    assert _close(tm.q_sca(), cs["q_sca"], tol=1e-9)
    print(f"    T-матрица q_sca={tm.q_sca():.6f} — согласовано")


def test_sphere_as_scatterer_in_cluster():
    print("\n[3] SphereSolver.as_scatterer работает в Cluster (GMM):")
    s1 = gt.SphereSolver(0.5, [2.25], position=(-1, 0, 0)).as_scatterer()
    s2 = gt.SphereSolver(0.5, [2.25], position=(+1, 0, 0)).as_scatterer()
    cl = gt.Cluster([s1, s2])
    cs_cluster = cl.cross_sections(k=3.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=5)
    cs_facade = gt.solve_cluster([s1, s2], 3.0, (0, 0, 1), (1, 0, 0), 5)
    for key in ("c_sca", "c_ext", "c_abs"):
        assert _close(cs_cluster[key], cs_facade[key], tol=1e-9), key
        assert np.isfinite(cs_cluster[key])
    print(f"    кластер C_sca={cs_cluster['c_sca']:.6f} (== facade) — ок")


def test_cluster_solve_shapes():
    print("\n[4] Cluster.solve возвращает согласованные (a, c, d):")
    s1 = gt.SphereSolver(0.5, [2.25], position=(0, 0, +1)).as_scatterer()
    s2 = gt.SphereSolver(0.5, [2.25], position=(0, 0, -1)).as_scatterer()
    a, c, d = gt.Cluster([s1, s2]).solve(k=3.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=4)
    ra, rc, rd = gmm.solve_cluster([s1, s2], 3.0, (0, 0, 1), (1, 0, 0), 4)
    blk = 2 * len(vswf.mode_list(4))   # 2K блоков (M,N), K = nmax(nmax+2)
    assert a.shape == ra.shape == (2, blk)
    assert np.allclose(a, ra) and np.allclose(c, rc) and np.allclose(d, rd)
    print(f"    формы a={a.shape}, совпадает с gmm.solve_cluster — ок")


def test_ellipsoid_solver_matches_module():
    print("\n[5] EllipsoidSolver == ellipsoid + solve_ellipsoid:")
    a, b, c, eps, k = 2.0, 1.0, 1.0, 3.0 + 0.1j, 0.5
    es = gt.EllipsoidSolver(a, b, c, eps)
    # per-axis
    al = ellipsoid.polarizability_homogeneous(a, b, c, eps, 1.0)
    ref0 = ellipsoid.rayleigh_cross_sections(al[0], k)
    got0 = es.cross_sections(k, axis=0)
    assert _close(got0["c_ext"], ref0["c_ext"]), "per-axis"
    # orientation average + facade
    refavg = ellipsoid.orientation_average_cross_sections(a, b, c, eps, k, 1.0)
    gotavg = es.cross_sections(k)
    facade = gt.solve_ellipsoid(a, b, c, eps, k)
    for key in ("c_sca", "c_abs", "c_ext"):
        assert _close(gotavg[key], refavg[key]), key
        assert _close(facade[key], refavg[key]), key
    assert np.allclose(es.polarizability(), al)
    print(f"    ⟨C_ext⟩={gotavg['c_ext']:.6f} (== module, == facade) — ок")


def test_spheroid_solver_matches_module():
    print("\n[6] SpheroidSolver == spheroid + solve_spheroid:")
    a_eq, c_ax, eps, k = 1.0, 2.0, 3.0 + 0.1j, 0.5
    ss = gt.SpheroidSolver(a_eq, c_ax, eps)
    assert ss.is_prolate is True
    L_eq_c, L_ax_c = spheroid.depolarization_closed(a_eq, c_ax)
    L = ss.depolarization_factors(closed=True)
    assert _close(L[0], L_eq_c) and _close(L[2], L_ax_c)
    refavg = ellipsoid.orientation_average_cross_sections(a_eq, a_eq, c_ax, eps, k, 1.0)
    facade = gt.solve_spheroid(a_eq, c_ax, eps, k)
    got = ss.cross_sections(k)
    for key in ("c_sca", "c_abs", "c_ext"):
        assert _close(got[key], refavg[key]), key
        assert _close(facade[key], refavg[key]), key
    print(f"    L=({L[0]:.4f},{L[1]:.4f},{L[2]:.4f}) ⟨C_ext⟩={got['c_ext']:.6f} — ок")


def test_cylinder_solver_matches_module():
    print("\n[7] CylinderSolver == cylinder + solve_cylinder:")
    radius, eps, k = 1.0, 4.0, 2.0
    cs = gt.CylinderSolver(radius, eps)
    assert _close(complex(cs.m).real, 2.0) and _close(complex(cs.m).imag, 0.0)
    ref = cylinder.cross_sections_infinite(2.0, k * radius, mode="TM")
    got = cs.cross_sections(k, mode="TM")
    facade = gt.solve_cylinder(radius, eps, k, mode="TM")
    for key in ("q_sca", "q_ext", "q_abs"):
        assert _close(got[key], ref[key]), key
        assert _close(facade[key], ref[key]), key
    print(f"    Q_sca={got['q_sca']:.6f} (== module, == facade) — ок")


def test_cone_solver_decompose_and_cluster():
    print("\n[8] ConeSolver.decompose -> Cluster (строгая аналитика через сферы):")
    cone_s = gt.ConeSolver(apex=[0, 0, 0], axis=[0, 0, 1],
                           half_angle=math.radians(35), height=8.0, eps=[2.25])
    scat, centers, radius = cone_s.decompose(spacing=2.0)
    assert len(scat) >= 1 and radius > 0
    cs = gt.Cluster(scat).cross_sections(k=1.0, khat=(1, 0, 0), pol=(0, 0, 1), nmax=3)
    assert all(np.isfinite(v) for v in cs.values())
    print(f"    сфер={len(scat)} radius={radius:.3f}; C_ext={cs['c_ext']:.5f} — ок")


def test_notimplemented_honesty():
    print("\n[9] полноволновые ветви → честный NotImplementedError:")
    raised = 0
    for fn in (lambda: gt.SpheroidSolver(1, 2, 2.0).full_wave(),
               lambda: gt.CylinderSolver(1, 2.0).finite(),
               lambda: gt.ConeSolver([0, 0, 0], [0, 0, 1], 0.5, 1.0, [2.0]).full_wave()):
        try:
            fn()
        except NotImplementedError:
            raised += 1
    assert raised == 3, "все три ветви должны поднимать NotImplementedError"
    print("    спфероид/конечный цилиндр/конус — 3/3 NotImplementedError — ок")


if __name__ == "__main__":
    for fn in (test_sphere_solver_matches_core, test_sphere_tmatrix_consistency,
               test_sphere_as_scatterer_in_cluster, test_cluster_solve_shapes,
               test_ellipsoid_solver_matches_module, test_spheroid_solver_matches_module,
               test_cylinder_solver_matches_module, test_cone_solver_decompose_and_cluster,
               test_notimplemented_honesty):
        fn()
    print("\nOK")
