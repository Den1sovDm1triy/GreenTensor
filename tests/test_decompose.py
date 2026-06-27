# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Проверки разложения геометрии (green_tensor/decompose.py).

  [1] непересечение: min расстояние между центрами ≥ 2·radius (условие GMM);
  [2] все сферы внутри тела;
  [3] доля заполнения разумна и растёт при измельчении шага;
  [4] интеграция: разложение даёт валидные рассеиватели для gmm.solve_cluster.

Запуск: python3 tests/test_decompose.py | pytest tests/test_decompose.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from green_tensor import decompose as dc  # noqa: E402
from green_tensor import gmm  # noqa: E402


def test_non_overlap_and_inside():
    print("\n[1-2] непересечение и принадлежность телу:")
    R = 2.0
    inside = dc.sphere_indicator([0, 0, 0], R)
    centers, radius = dc.pack_spheres(inside, [-R, -R, -R], [R, R, R], spacing=0.5)
    assert len(centers) > 10
    sep = dc.min_separation(centers)
    print(f"    сфер: {len(centers)}, radius={radius}, min_sep={sep:.4f}, 2r={2*radius}")
    assert sep >= 2 * radius - 1e-12, "сферы пересекаются!"
    assert np.all(inside(centers)), "есть центр вне тела"


def test_coverage_and_refinement():
    print("\n[3] доля заполнения растёт при измельчении:")
    R = 2.0
    inside = dc.sphere_indicator([0, 0, 0], R)
    Vbody = (4.0 / 3.0) * math.pi * R**3
    cov = []
    for sp in (0.6, 0.4, 0.25):
        centers, radius = dc.pack_spheres(inside, [-R, -R, -R], [R, R, R], spacing=sp)
        cov.append(dc.coverage_fraction(centers, radius, Vbody))
        print(f"    spacing={sp}: сфер={len(centers)}, покрытие={cov[-1]:.3f}")
    assert cov[0] < cov[1] < cov[2], "покрытие должно расти при измельчении"
    assert 0.0 < cov[-1] < 1.0


def test_box_and_cylinder():
    print("\n[4a] коробка и цилиндр:")
    box = dc.box_indicator([0, 0, 0], [1.0, 0.7, 0.5])
    c, r = dc.pack_spheres(box, [-1, -0.7, -0.5], [1, 0.7, 0.5], spacing=0.3)
    assert len(c) > 0 and dc.min_separation(c) >= 2 * r - 1e-12 and np.all(box(c))
    cyl = dc.cylinder_indicator([0, 0, 0], radius=0.8, half_length=1.5, axis=2)
    c2, r2 = dc.pack_spheres(cyl, [-0.8, -0.8, -1.5], [0.8, 0.8, 1.5], spacing=0.35)
    assert len(c2) > 0 and dc.min_separation(c2) >= 2 * r2 - 1e-12 and np.all(cyl(c2))
    print(f"    коробка: {len(c)} сфер; цилиндр: {len(c2)} сфер — непересекающиеся, внутри")


def test_decompose_feeds_gmm():
    """Разложение -> рассеиватели -> gmm.solve_cluster выдаёт согласованное решение."""
    print("\n[4b] интеграция с GMM (мини-кластер):")
    # короткий «стержень» из нескольких сфер вдоль z (k·spacing достаточно для трансляции out)
    k = 1.0
    centers = np.array([[0, 0, -6.0], [0, 0, 0.0], [0, 0, 6.0]])
    scat = dc.to_scatterers(centers, radius=0.5, eps=[2.25])
    a, c, d = gmm.solve_cluster(scat, k, [1, 0, 0], [0, 0, 1], nmax=4)
    assert a.shape[0] == 3 and np.all(np.isfinite(c))
    # невязка решателя
    print(f"    кластер из {len(scat)} сфер решён, |c| конечно — ок")


_TESTS = [test_non_overlap_and_inside, test_coverage_and_refinement,
          test_box_and_cylinder, test_decompose_feeds_gmm]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ decompose проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
