# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Проверки decompose-fallback конуса (green_tensor/decompose.py: cone_indicator,
decompose_cone). Строгий полноволновой конус — EBCM-примитив green_tensor.Cone
(см. tests/test_ebcm.py::test_ebcm_cone_energy).

  [1] индикатор конуса: точки внутри/снаружи классифицируются верно;
  [2] разложение конуса в сферы: непересечение + все внутри + интеграция с GMM.

Запуск: python3 tests/test_cone.py | pytest tests/test_cone.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from green_tensor import decompose as cone  # noqa: E402  (cone_indicator/decompose_cone теперь в decompose)
from green_tensor import decompose as dc  # noqa: E402
from green_tensor import gmm  # noqa: E402


def test_cone_indicator():
    print("\n[1] индикатор конуса:")
    f = cone.cone_indicator([0, 0, 0], [0, 0, 1], half_angle=math.radians(30), height=2.0)
    pts = np.array([
        [0.0, 0.0, 1.0],     # на оси внутри
        [0.5, 0.0, 1.0],     # rad=0.5 ≤ 1*tan30≈0.577 — внутри
        [0.7, 0.0, 1.0],     # rad=0.7 > 0.577 — снаружи
        [0.0, 0.0, 2.5],     # за высотой — снаружи
        [0.0, 0.0, -0.1],    # за вершиной — снаружи
    ])
    res = f(pts)
    print(f"    {res.tolist()}")
    assert res.tolist() == [True, True, False, False, False]


def test_decompose_cone_and_gmm():
    print("\n[2] разложение конуса в сферы + GMM:")
    scat, centers, radius = cone.decompose_cone(
        apex=[0, 0, 0], axis=[0, 0, 1], half_angle=math.radians(35),
        height=8.0, spacing=2.0, eps=[2.25])
    assert len(scat) >= 2, f"слишком мало сфер: {len(scat)}"
    sep = dc.min_separation(centers)
    inside = cone.cone_indicator([0, 0, 0], [0, 0, 1], math.radians(35), 8.0)
    print(f"    сфер={len(scat)} radius={radius:.3f} min_sep={sep:.3f} 2r={2*radius:.3f}")
    assert sep >= 2 * radius - 1e-12 and np.all(inside(centers))
    a, c, d = gmm.solve_cluster(scat, k=1.0, khat=[1, 0, 0], pol=[0, 0, 1], nmax=3)
    assert np.all(np.isfinite(c))
    print(f"    GMM на конус-кластере решён, |c| конечно — ок")


_TESTS = [test_cone_indicator, test_decompose_cone_and_gmm]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ cone проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
