# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Сложная геометрия из примитивов (CSG-объединение) через когерентный GMM.

Тело собирается как набор РАЗМЕЩЁННЫХ аналитических примитивов (сфера/сфероид/конус/
конечный цилиндр), каждый даёт сферическую T-матрицу (сфера — Ми, остальные — EBCM),
а GMM решает связанную систему. Проверяем:
  [1] редукция: кластер из ОДНОГО примитива == одиночный примитив (бит-в-бит);
  [2] энергобаланс составного тела без потерь (C_abs≈0) для близких примитивов;
  [3] расцепление: при большом разносе коэффициенты тела → одиночные (связь → 0).

Ось симметрии примитивов ∥ z (произвольная ориентация — через вращение Вигнера-D,
ещё не реализовано); положение — произвольное (трансляция GMM).

Запуск: python3 tests/test_composite.py | pytest tests/test_composite.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from green_tensor import gmm  # noqa: E402
from green_tensor.scatterer import LayeredSphere, Spheroid, Cone, FiniteCylinder  # noqa: E402

K = 2.0
NMAX = 7
KHAT = (1.0, 0.0, 0.0)
POL = (0.0, 0.0, 1.0)


def _cluster_cs(scatterers):
    _, c, _ = gmm.solve_cluster(scatterers, K, KHAT, POL, NMAX)
    return gmm.cross_sections(scatterers, c, K, KHAT, POL)


def test_single_primitive_reduction():
    """Кластер из одного примитива == одиночное решение (gmm.isolated_scattered)."""
    print("\n[1] редукция: Cluster([X]) == одиночный X:")
    prims = [("spheroid", Spheroid((0, 0, 0), 0.4, 0.8, 2.25)),
             ("cone", Cone((0, 0, 0), 0.4, 0.8, 2.25)),
             ("cylinder", FiniteCylinder((0, 0, 0), 0.3, 0.8, 2.25))]
    for name, p in prims:
        c_iso = gmm.isolated_scattered(p, K, KHAT, POL, NMAX)
        cs_iso = gmm.cross_sections([p], [c_iso], K, KHAT, POL)
        cs_cl = _cluster_cs([p])
        d = abs(cs_cl["c_sca"] - cs_iso["c_sca"]) / cs_iso["c_sca"]
        print(f"    {name}: |ΔC_sca|/C_sca = {d:.1e}")
        assert d < 1e-12, f"{name}: кластер-из-одного должен совпасть с одиночным"


def test_composite_energy_conservation():
    """Составное тело без потерь: C_abs = C_ext − C_sca ≈ 0 (близкие примитивы)."""
    print("\n[2] энергобаланс составного тела (lossless):")
    cases = [("sphere+spheroid", [LayeredSphere((0, 0, 1.2), 0.4, [2.25]),
                                   Spheroid((0, 0, -1.0), 0.35, 0.7, 2.25)], 5e-3),
             ("sphere+cone", [LayeredSphere((0, 0, 1.3), 0.4, [2.25]),
                              Cone((0, 0, -0.8), 0.4, 0.8, 2.25)], 3e-2),
             ("two spheroids", [Spheroid((0, 0, 1.0), 0.35, 0.6, 2.25),
                                Spheroid((0, 0, -1.0), 0.35, 0.6, 2.25)], 5e-3)]
    for name, comb, tol in cases:
        cs = _cluster_cs(comb)
        bal = abs(cs["c_abs"]) / cs["c_sca"]
        print(f"    {name}: C_sca={cs['c_sca']:.4f} |C_abs|/C_sca={bal:.1e} (tol={tol:.0e})")
        assert bal < tol, f"{name}: энергобаланс составного тела нарушен ({bal:.1e})"


def test_composite_decoupling():
    """Расцепление: при росте разноса коэффициенты примитива → одиночные (связь → 0)."""
    print("\n[3] расцепление составного тела (коэффициенты → одиночные):")
    errs = []
    for sep in (4.0, 40.0, 400.0):
        s1 = Spheroid((0, 0, +sep / 2), 0.4, 0.8, 2.25)
        s2 = LayeredSphere((0, 0, -sep / 2), 0.4, [2.25])
        _, c, _ = gmm.solve_cluster([s1, s2], K, KHAT, POL, NMAX)
        c_iso = gmm.isolated_scattered(s1, K, KHAT, POL, NMAX)
        e = np.max(np.abs(c[0] - c_iso)) / np.max(np.abs(c_iso))
        errs.append(e)
        print(f"    k·sep={K*sep:5.0f}: ||c − c_iso||/||c_iso|| = {e:.2e}")
    assert errs[0] > errs[1] > errs[2], "ошибка связи должна убывать с разносом"
    assert errs[2] < 1e-4, f"при большом разносе должно сойтись к одиночному ({errs[2]:.1e})"


_TESTS = [test_single_primitive_reduction, test_composite_energy_conservation,
          test_composite_decoupling]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ composite проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
