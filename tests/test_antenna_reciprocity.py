# SPDX-License-Identifier: MIT
"""Антенный режим через теорему взаимности: H-поле и ДН линзовой антенны.

Проверяется:
  1. аналитическое H совпадает с ∇×E/(jkZ0) численного дифференцирования
     аналитического E — в обоих направлениях падения, рассеянное и полное;
  2. в дальней зоне H = r̂ × E / Z0 (соотношение исходящей волны);
  3. комбинация Гюйгенса в свободном пространстве даёт ТОЧНО кардиоиду
     (1+cosΘ)/2 с максимумом внутрь тела — этим закреплены знак и константа;
  4. ДН линзы Люнеберга k0a=15 совпадает с эталоном Ansys HFSS.
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

from green_tensor.nearfield import (  # noqa: E402
    antenna_pattern_sphere, external_field_sphere, magnetic_field_sphere)

Z0 = 376.730313668
F0, C0 = 10e9, 299792458.0
K0 = 2.0 * math.pi * F0 / C0
LAM = C0 / F0
A_LENS = 71.5702e-3
EPS = [1.96, 1.84, 1.64, 1.36, 1.00]
ANORM = [0.2, 0.4, 0.6, 0.8, 1.0]
GAP = 0.074948115e-3


def _curl_H(pts, khat, total):
    """H = ∇×E/(j k Z0) центральными разностями по аналитическому E."""
    kw = dict(problem="diffraction", a_norm=ANORM, khat=khat, e0=1.0,
              mask_inside=False, total=total)
    p = np.asarray(pts, float)
    step = LAM / 4000.0
    d = np.zeros((3, 3, len(p)), dtype=complex)
    for i in range(3):
        sh = np.zeros(3); sh[i] = step
        d[i] = ((external_field_sphere(A_LENS, EPS, K0, p + sh, **kw)
                 - external_field_sphere(A_LENS, EPS, K0, p - sh, **kw)) / (2 * step)).T
    curl = np.stack([d[1][2] - d[2][1], d[2][0] - d[0][2], d[0][1] - d[1][0]], axis=1)
    return curl / (1j * K0 * Z0)


def test_magnetic_field_matches_curl():
    rng = np.random.default_rng(0)
    pts = []
    for _ in range(6):
        r = A_LENS * (1.02 + 2.0 * rng.random())
        t, f = rng.uniform(0.2, 2.9), rng.uniform(0.0, 2 * math.pi)
        pts.append([r * math.sin(t) * math.cos(f), r * math.sin(t) * math.sin(f),
                    r * math.cos(t)])
    pts = np.array(pts)
    for khat in ((0.0, 0.0, -1.0), (0.0, 0.0, 1.0)):
        for total in (False, True):
            Ha = magnetic_field_sphere(A_LENS, EPS, K0, pts, a_norm=ANORM, khat=khat,
                                       e0=1.0, total=total, mask_inside=False)
            Hc = _curl_H(pts, khat, total)
            rel = np.max(np.linalg.norm(Ha - Hc, axis=1) / np.linalg.norm(Hc, axis=1))
            assert rel < 1e-5, f"khat={khat} total={total}: {rel:.2e}"


def test_far_zone_outgoing_relation():
    th = np.deg2rad(np.array([30.0, 75.0, 120.0, 165.0]))
    R = 60.0 * LAM
    pts = np.stack([R * np.sin(th), np.zeros_like(th), R * np.cos(th)], axis=1)
    kw = dict(a_norm=ANORM, khat=(0.0, 0.0, -1.0), e0=1.0, mask_inside=False)
    E = external_field_sphere(A_LENS, EPS, K0, pts, **kw)
    H = magnetic_field_sphere(A_LENS, EPS, K0, pts, **kw)
    rh = pts / np.linalg.norm(pts, axis=1)[:, None]
    rel = np.linalg.norm(H - np.cross(rh, E) / Z0, axis=1) / np.linalg.norm(E / Z0, axis=1)
    assert np.max(rel) < 2e-3, f"{np.max(rel):.2e}"


def test_huygens_free_space_is_cardioid():
    """Свободное пространство: комбинация Гюйгенса == 2·(1+cosΘ)/2, точно."""
    th = np.deg2rad(np.arange(1.0, 360.0, 1.0))
    rd = A_LENS + GAP
    pts = np.stack([rd * np.sin(th), np.zeros_like(th), rd * np.cos(th)], axis=1)
    kw = dict(a_norm=ANORM, khat=(0.0, 0.0, -1.0), e0=1.0, mask_inside=False)
    Ei = (external_field_sphere(A_LENS, EPS, K0, pts, total=True, **kw)
          - external_field_sphere(A_LENS, EPS, K0, pts, total=False, **kw))
    Hi = (magnetic_field_sphere(A_LENS, EPS, K0, pts, total=True, **kw)
          - magnetic_field_sphere(A_LENS, EPS, K0, pts, total=False, **kw))
    pe = np.stack([np.cos(th), np.zeros_like(th), -np.sin(th)], axis=1)
    pm = np.tile(np.array([0.0, 1.0, 0.0]), (len(th), 1))
    w = np.abs(np.sum(Ei * pe, axis=1) + Z0 * np.sum(Hi * pm, axis=1))
    assert abs(np.max(w) - 2.0) < 1e-9, f"пик {np.max(w)!r} != 2"
    card = 2.0 * (1.0 + np.cos(math.pi - th)) / 2.0     # максимум внутрь тела
    assert np.max(np.abs(w - card)) < 1e-9, f"{np.max(np.abs(w - card)):.2e}"


def test_luneburg_antenna_vs_hfss():
    """Линза Люнеберга k0a=15: ДН по взаимности против эталона Ansys HFSS."""
    ref = os.path.join(_ROOT, "examples",
                       "Example 7 - Absolute fields and lens antenna",
                       "reference", "csv", "hfss_far_luneburg_k15_antenna.csv")
    if not os.path.exists(ref):
        return                                       # эталон не входит в пакет
    h = np.loadtxt(ref, delimiter=",", skiprows=1)
    hth, hv = h[:, 0], h[:, 1] / 10.0
    th = np.deg2rad(np.arange(0.5, 360.0, 0.5))
    out = antenna_pattern_sphere(A_LENS, EPS, K0, th, current=1.0, length=LAM / 50.0,
                                 gap=GAP, r=10.0, plane="E", a_norm=ANORM)
    v, thd = out["E"], np.degrees(th)
    # абсолютный уровень
    d_peak = 20 * math.log10(np.max(v) / np.max(hv))
    assert abs(d_peak) < 1.0, f"Δпик = {d_peak:+.2f} дБ"
    assert abs(thd[int(np.argmax(v))] - 180.0) < 1.0
    # форма по всей диаграмме
    hi = np.interp(thd % 360, hth % 360, hv, period=360)
    m = v > 0.01 * np.max(v)
    db = lambda x: 20 * np.log10(np.maximum(x, 1e-300))
    med = float(np.nanmedian(np.abs(db(hi[m] / np.max(hi)) - db(v[m] / np.max(v)))))
    assert med < 1.0, f"медиана формы {med:.2f} дБ"


if __name__ == "__main__":
    for fn in (test_magnetic_field_matches_curl, test_far_zone_outgoing_relation,
               test_huygens_free_space_is_cardioid, test_luneburg_antenna_vs_hfss):
        fn(); print(f"  ok   {fn.__name__}")
    print("ВСЕ ТЕСТЫ ПРОЙДЕНЫ")
