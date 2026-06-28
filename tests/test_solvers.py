# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Проверки публичного API green_tensor/solvers.py.

Каждый класс-решатель и функция-обёртка ``solve_*`` НЕ содержат своей математики —
тесты подтверждают, что они в точности воспроизводят уже проверенные нижележащие
функции (никакого расхождения), и что полноволновые ветви честно поднимают
NotImplementedError.
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

import green_tensor as gt  # noqa: E402
from analytic_mie import _q_ext, _q_sca, mie_ab_eps_mu  # noqa: E402
from green_tensor import cylinder, gmm, sphere_core, tmatrix, vswf  # noqa: E402


def _close(a, b, tol=1e-12):
    return abs(a - b) <= tol * (1.0 + abs(b))


def test_sphere_solver_matches_core():
    print("\n[1] SphereSolver == canonical sphere facade == solve_sphere:")
    radius, eps, k = 1.0, [2.25], 3.0
    sp = gt.SphereSolver(radius, eps)
    cs = sp.cross_sections(k)
    ref = sphere_core.MieSphere(k0=k * radius, a=[1.0], eps=eps).cross_sections()
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


def test_notimplemented_honesty():
    print("\n[6] нереализованная ветвь → честный NotImplementedError:")
    # Конечный цилиндр несепарабелен и не имеет точного аналитического решения в этом
    # базисе; заглушка честно поднимает ошибку (точна только бесконечная аналитика).
    raised = False
    try:
        gt.CylinderSolver(1, 2.0).finite()
    except NotImplementedError:
        raised = True
    assert raised, "cylinder.finite должен поднимать NotImplementedError"
    print("    cylinder.finite → NotImplementedError; точна только бесконечная аналитика — ок")


def test_sphere_antenna_problem():
    print("\n[7] sphere antenna problem (источник на поверхности) — канон 01_sphere:")
    sp = gt.SphereSolver(1.0, [2.25])
    pd = sp.pattern(3.0, problem="diffraction")["E"]
    pa = sp.pattern(3.0, problem="antenna")["E"]
    assert np.all(np.isfinite(pd)) and np.all(np.isfinite(pa)), "antenna pattern must be finite"
    assert not np.allclose(pd, pa), "antenna must differ from diffraction"
    p0 = sp.pattern(3.0)["E"]                       # default branch
    assert np.allclose(p0, pd), "default problem must be diffraction (backward compat)"
    print(f"    antenna max|E|={np.nanmax(pa):.3f} ≠ diffraction {np.nanmax(pd):.3f}; default==diffraction — ок")


if __name__ == "__main__":
    for fn in (test_sphere_solver_matches_core, test_sphere_tmatrix_consistency,
               test_sphere_solver_magnetodielectric_matches_closed_mie,
               test_sphere_as_scatterer_in_cluster, test_cluster_solve_shapes,
               test_cylinder_solver_matches_module, test_notimplemented_honesty,
               test_sphere_antenna_problem):
        fn()
    print("\nOK")
