# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

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
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

import green_tensor as gt  # noqa: E402
from analytic_mie import _q_ext, _q_sca, mie_ab_eps_mu  # noqa: E402
from green_tensor import cylinder, ellipsoid, gmm, mie_core, tmatrix, vswf  # noqa: E402


def _close(a, b, tol=1e-12):
    return abs(a - b) <= tol * (1.0 + abs(b))


def test_sphere_solver_matches_core():
    print("\n[1] SphereSolver == canonical sphere facade == solve_sphere:")
    radius, eps, k = 1.0, [2.25], 3.0
    sp = gt.SphereSolver(radius, eps)
    cs = sp.cross_sections(k)
    ref = mie_core.MieSphere(k0=k * radius, a=[1.0], eps=eps).cross_sections()
    facade = gt.solve_sphere(radius, eps, k)
    for key in ("q_sca", "q_ext", "q_abs", "q_back"):
        assert _close(cs[key], ref[key]), key
        assert _close(facade[key], ref[key]), key
    print(f"    Q_sca={cs['q_sca']:.6f} (== canonical, == facade) — ок")


def test_sphere_tmatrix_consistency():
    print("\n[2] SphereSolver.t_matrix.q_sca == cross_sections:")
    sp = gt.SphereSolver(1.0, [2.25])
    cs = sp.cross_sections(3.0)
    tm = sp.t_matrix(3.0)
    assert isinstance(tm, tmatrix.DiagonalTMatrix)
    assert _close(tm.q_sca(), cs["q_sca"], tol=1e-9)
    print(f"    T-матрица q_sca={tm.q_sca():.6f} — согласовано")


def test_sphere_solver_magnetodielectric_matches_closed_mie():
    print("\n[2b] SphereSolver магнитодиэлектрик eps,mu == закрытая формула Ми:")
    eps, mu, x, toch = 1.0, 4.0, 0.5, 10
    sp = gt.SphereSolver(1.0, [eps], miy=[mu])
    got = sp.cross_sections(x, toch=toch)
    facade = gt.solve_sphere(1.0, [eps], x, miy=[mu], toch=toch)
    tm = sp.t_matrix(x, toch=toch)
    n, a, b = mie_ab_eps_mu(eps, mu, x, toch)
    ref_sca, ref_ext = _q_sca(n, a, b, x), _q_ext(n, a, b, x)
    assert _close(got["q_sca"], ref_sca, tol=1e-10)
    assert _close(got["q_ext"], ref_ext, tol=1e-10)
    assert _close(facade["q_sca"], ref_sca, tol=1e-10)
    assert _close(tm.q_sca(), ref_sca, tol=1e-10)
    print(f"    Q_sca={got['q_sca']:.6e}, Q_ext={got['q_ext']:.6e} — ок")


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


def test_spheroid_primitive_ebcm_in_cluster():
    print("\n[6] gt.Spheroid (строгий EBCM) в Cluster == gmm одиночным телом:")
    a_eq, c_ax, eps, k = 0.4, 0.8, 2.25, 2.0
    sp = gt.Spheroid((0, 0, 0), a_eq, c_ax, eps)
    khat, pol, nmax = (0, 0, 1), (1, 0, 0), 6
    got = gt.Cluster([sp]).cross_sections(k=k, khat=khat, pol=pol, nmax=nmax)
    c = gmm.isolated_scattered(sp, k, khat, pol, nmax)
    ref = gmm.cross_sections([sp], np.array([c]), k, khat, pol)
    for key in ("c_sca", "c_ext", "c_abs"):
        assert _close(got[key], ref[key], tol=1e-9), key
        assert np.isfinite(got[key])
    bal = abs(got["c_abs"]) / got["c_sca"]
    assert bal < 5e-3, f"диэлектрический сфероид: энергобаланс нарушен ({bal:.1e})"
    print(f"    C_sca={got['c_sca']:.6f} |C_abs|/C_sca={bal:.1e} (lossless) — ок")


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
    print("\n[9] нереализованная ветвь → честный NotImplementedError:")
    # Конечный цилиндр в 2D-модуле cylinder.py намеренно не реализован (мод-матчинг);
    # строгий путь — полноволновой EBCM gt.FiniteCylinder. Заглушка честно поднимает ошибку.
    raised = False
    try:
        gt.CylinderSolver(1, 2.0).finite()
    except NotImplementedError:
        raised = True
    assert raised, "cylinder.finite (2D-модуль) должен поднимать NotImplementedError"
    print("    cylinder.finite (2D) → NotImplementedError; строгий конечный цилиндр = gt.FiniteCylinder — ок")


if __name__ == "__main__":
    for fn in (test_sphere_solver_matches_core, test_sphere_tmatrix_consistency,
               test_sphere_solver_magnetodielectric_matches_closed_mie,
               test_sphere_as_scatterer_in_cluster, test_cluster_solve_shapes,
               test_ellipsoid_solver_matches_module, test_spheroid_primitive_ebcm_in_cluster,
               test_cylinder_solver_matches_module, test_cone_solver_decompose_and_cluster,
               test_notimplemented_honesty):
        fn()
    print("\nOK")
