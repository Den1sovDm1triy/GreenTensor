# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Вращение Вигнера-D: произвольная ориентация рассеивателя (поворот T-матрицы).

Тело строится в собственной системе (ось ∥ z); ориентация задаётся углами Эйлера
(α,β,γ), z-y-z, активный поворот R=Rz(α)Ry(β)Rz(γ); T_global = D(α,β,γ)·T·D^†.
Проверяем:
  [1] малая d-матрица Вигнера: замкнутая форма d^1, ортогональность, композиция;
  [2] T сферы инвариантна к любому повороту; T осесимметричного тела — к повороту
      вокруг своей оси (β=0);
  [3] КОВАРИАНТНОСТЬ: сечение C_sca инвариантно при одновременном повороте тела,
      падения и поляризации на один R — для произвольных (α,β,γ);
  [4] ориентированный примитив в GMM (наклонный сфероид): энергобаланс сохраняется.

Запуск: python3 tests/test_rotation.py | pytest tests/test_rotation.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from green_tensor import ebcm, gmm, vswf  # noqa: E402
from green_tensor.scatterer import Spheroid, LayeredSphere  # noqa: E402

K = 2.0
NMAX = 7


def _Rzyz(a, b, g):
    def Rz(t):
        return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])

    def Ry(t):
        return np.array([[np.cos(t), 0, np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])
    return Rz(a) @ Ry(b) @ Rz(g)


class _TW:
    def __init__(self, T):
        self._T = T
        self.position = np.zeros(3)

    def t_matrix(self, k, nmax):
        return self._T


def _csca(T, khat, pol):
    s = _TW(T)
    c = gmm.isolated_scattered(s, K, khat, pol, NMAX)
    return gmm.cross_sections([s], [c], K, khat, pol)["c_sca"]


def test_wigner_d_closed_form():
    print("\n[1] малая d-матрица Вигнера:")
    b = 0.5
    c, s = np.cos(b), np.sin(b)
    d = vswf.wigner_d(1, b)            # индексы [m'+1, m+1]
    # стандарт: d^1_{1,0}=−s/√2 → d[2,1]; d^1_{0,0}=c → d[1,1]; d^1_{1,1}=(1+c)/2 → d[2,2]
    ok = (abs(d[2, 1] - (-s / np.sqrt(2))) < 1e-12 and abs(d[1, 1] - c) < 1e-12
          and abs(d[2, 2] - (1 + c) / 2) < 1e-12 and abs(d[1, 2] - (s / np.sqrt(2))) < 1e-12)
    print(f"    d^1 vs замкнутая форма: {'OK' if ok else 'FAIL'}")
    assert ok, "d^1 не совпала с замкнутой формой"
    d4 = vswf.wigner_d(4, 0.7)
    assert np.allclose(d4 @ d4.T, np.eye(9)), "d не ортогональна"
    assert np.allclose(vswf.wigner_d(3, 0.3) @ vswf.wigner_d(3, 0.4), vswf.wigner_d(3, 0.7)), \
        "композиция d(β1)d(β2)≠d(β1+β2)"
    print("    ортогональность + композиция: OK")


def test_sphere_and_axisym_invariance():
    print("\n[2] инвариантность T при повороте:")
    Tsph = ebcm.tmatrix_axisym(ebcm.sphere_curve(0.5), K, 2.25, 1.0, NMAX)
    d = np.max(np.abs(vswf.rotate_tmatrix(Tsph, NMAX, 0.5, 0.9, 0.7) - Tsph))
    print(f"    сфера инвариантна к любому повороту: max|ΔT|={d:.1e}")
    assert d < 1e-12, "T сферы должна быть инвариантна к повороту"
    Tspd = ebcm.tmatrix_axisym(ebcm.spheroid_curve(0.4, 0.8), K, 2.25, 1.0, NMAX)
    d2 = np.max(np.abs(vswf.rotate_tmatrix(Tspd, NMAX, 1.1, 0.0, 0.0) - Tspd))
    print(f"    сфероид инвариантен к повороту вокруг своей оси (β=0): max|ΔT|={d2:.1e}")
    assert d2 < 1e-12, "T осесимметричного тела должна быть инвариантна к Rz"


def test_rotation_covariance():
    print("\n[3] ковариантность сечения (поворот тела+падения+поляризации):")
    T = ebcm.tmatrix_axisym(ebcm.spheroid_curve(0.4, 0.8), K, 2.25, 1.0, NMAX)
    th0, ph0 = np.radians(40.0), np.radians(25.0)
    khat0 = np.array([np.sin(th0) * np.cos(ph0), np.sin(th0) * np.sin(ph0), np.cos(th0)])
    pol0 = np.cross(khat0, [0, 0, 1.0])
    pol0 = pol0 / np.linalg.norm(pol0)
    C_A = _csca(T, khat0, pol0)
    for a, b, g in [(0.0, 0.6, 0.0), (0.5, 0.9, 0.7), (1.2, 1.3, 2.0), (2.5, 0.4, 1.0)]:
        R = _Rzyz(a, b, g)
        Tr = vswf.rotate_tmatrix(T, NMAX, a, b, g)
        rel = abs(_csca(Tr, R @ khat0, R @ pol0) - C_A) / C_A
        print(f"    (α,β,γ)=({a},{b},{g}): rel={rel:.1e}")
        assert rel < 1e-5, f"ковариантность нарушена при (α,β,γ)=({a},{b},{g}): {rel:.1e}"


def test_oriented_scatterer_energy():
    print("\n[4] ориентированный примитив в GMM (наклонные сфероиды), энергобаланс:")
    khat, pol = (1.0, 0.0, 0.0), (0.0, 0.0, 1.0)
    body = [Spheroid((0, 0, 1.0), 0.35, 0.7, 2.25, euler=(0.0, 0.7, 0.0)),
            LayeredSphere((0, 0, -1.0), 0.4, [2.25])]
    _, c, _ = gmm.solve_cluster(body, K, khat, pol, NMAX)
    cs = gmm.cross_sections(body, c, K, khat, pol)
    bal = abs(cs["c_abs"]) / cs["c_sca"]
    print(f"    tilted spheroid + sphere: C_sca={cs['c_sca']:.4f} |C_abs|/C_sca={bal:.1e}")
    assert bal < 5e-3, f"энергобаланс с ориентированным телом нарушен ({bal:.1e})"
    # ориентированный сфероид t_matrix == rotate_tmatrix(неориентированный)
    sp_o = Spheroid((0, 0, 0), 0.4, 0.8, 2.25, euler=(0.3, 0.7, 1.1))
    sp_0 = Spheroid((0, 0, 0), 0.4, 0.8, 2.25)
    d = np.max(np.abs(sp_o.t_matrix(K, NMAX)
                      - vswf.rotate_tmatrix(sp_0.t_matrix(K, NMAX), NMAX, 0.3, 0.7, 1.1)))
    print(f"    Spheroid(euler).t_matrix == rotate_tmatrix: max|Δ|={d:.1e}")
    assert d < 1e-12


_TESTS = [test_wigner_d_closed_form, test_sphere_and_axisym_invariance,
          test_rotation_covariance, test_oriented_scatterer_energy]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ rotation проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
