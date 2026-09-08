# SPDX-License-Identifier: MIT

"""Верификация абсолютных амплитуд рассеяния против опубликованного эталона.

Источник эталона — стандартный тестовый набор MIEV0 У. Дж. Уискомба:
    W. J. Wiscombe, "Improved Mie scattering algorithms",
    Applied Optics 19(9), 1505–1509 (1980);
    NCAR Technical Note NCAR/TN-140+STR (тест-драйвер MVTstNew.f, зашитые эталоны).
Это универсальный бенчмарк для программ рассеяния Ми: таблица содержит АБСОЛЮТНЫЕ
комплексные амплитудные функции S_1(θ), S_2(θ) и эффективности Q_ext, Q_sca.

Проверяются оба звена: независимый ряд Ми (``analytic_mie``) и АБСОЛЮТНЫЕ амплитуды
нового ядра (``bistatic_rcs_sphere`` / ``MieSphere.amplitudes``). Совпадение — до
точности печати таблицы (7 значащих цифр, ~1e-7).

Конвенция: эталон Уискомба записывает поглощающий показатель как m = n − i·k; в
конвенции библиотеки (e^{-iωt}, Im ε > 0) та же среда — m = n + i·k, ε = m². Модули
|S_1|, |S_2| и эффективности от выбора конвенции не зависят.

Запуск:
    python3 tests/test_wiscombe_benchmark.py     # standalone
    pytest tests/test_wiscombe_benchmark.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from green_tensor import SphereSolver, bistatic_rcs_sphere  # noqa: E402
import analytic_mie as am  # noqa: E402

# Углы таблицы: θ = 0,30,…,180°. Значения (S1, S2) — verbatim из MVTstNew.f.
_ANGLES = np.array([0, 30, 60, 90, 120, 150, 180])

# --- Кейс 14 MIEV0: x=1.0, m=1.5−i (классический пример Bohren & Huffman) ---
CASE_14 = dict(
    name="Wiscombe 14: x=1.0, m=1.5-i (B&H)",
    x=1.0, m=1.5 + 1.0j, q_ext=2.336321, q_sca=6.634538e-01,
    S1=[5.840802e-01 + 1.905153e-01j, 5.657020e-01 + 1.871997e-01j,
        5.175251e-01 + 1.784426e-01j, 4.563396e-01 + 1.671665e-01j,
        4.002117e-01 + 1.566427e-01j, 3.621572e-01 + 1.493910e-01j,
        3.488438e-01 + 1.468286e-01j],
    S2=[5.840802e-01 + 1.905153e-01j, 5.001610e-01 + 1.456112e-01j,
        2.879639e-01 + 4.105398e-02j, 3.622847e-02 - 6.182646e-02j,
        -1.748750e-01 - 1.229586e-01j, -3.056823e-01 - 1.438460e-01j,
        -3.488438e-01 - 1.468286e-01j],
)

# --- Кейс 7 MIEV0: x=10, m=0.75 (непоглощающий, электрически крупнее) ---
CASE_7 = dict(
    name="Wiscombe 7: x=10, m=0.75",
    x=10.0, m=0.75 + 0.0j, q_ext=2.232265, q_sca=2.232265,
    S1=[5.580662e+01 - 9.758097e+00j, -7.672879e+00 + 1.087317e+01j,
        3.587894e+00 - 1.756177e+00j, -1.785905e+00 - 5.232828e-02j,
        1.537971e+00 - 8.329374e-02j, -4.140427e-01 + 1.876851e-01j,
        -1.078568e+00 - 3.608807e-02j],
    S2=[5.580662e+01 - 9.758097e+00j, -1.092923e+01 + 9.629667e+00j,
        3.427411e+00 + 8.082691e-02j, -5.148748e-01 - 7.027288e-01j,
        -6.908338e-01 + 2.152693e-01j, 5.247557e-01 - 1.923391e-01j,
        1.078568e+00 + 3.608807e-02j],
)

_CASES = [CASE_14, CASE_7]
_RTOL = 2e-6  # предел точности печати таблицы (7 знаков)


def _evaluate(case):
    """Вернуть (|S1|,|S2|) ссылки/Ми/ядра, (Qext,Qsca) ядра — для одного кейса."""
    x, eps = case["x"], case["m"] ** 2
    th = np.deg2rad(_ANGLES.astype(float))
    ref1 = np.abs(np.asarray(case["S1"]))
    ref2 = np.abs(np.asarray(case["S2"]))
    mie1, mie2 = am.mie_amplitudes(case["m"], x, th)
    solver = SphereSolver(1.0, [eps])          # R=1, k=x → k·R=x
    core1, core2 = solver.mie(x).amplitudes(th)
    cs = solver.cross_sections(k=x)
    return dict(ref1=ref1, ref2=ref2, mie1=np.abs(mie1), mie2=np.abs(mie2),
                core1=np.abs(core1), core2=np.abs(core2),
                q_ext=cs["q_ext"], q_sca=cs["q_sca"], q_back=cs["q_back"], x=x)


def test_amplitudes_match_wiscombe():
    """|S1|,|S2| ядра и ряда Ми совпадают с таблицей Уискомба (~1e-7)."""
    for case in _CASES:
        r = _evaluate(case)
        np.testing.assert_allclose(r["mie1"], r["ref1"], rtol=_RTOL, err_msg=case["name"] + " Mie |S1|")
        np.testing.assert_allclose(r["mie2"], r["ref2"], rtol=_RTOL, err_msg=case["name"] + " Mie |S2|")
        np.testing.assert_allclose(r["core1"], r["ref1"], rtol=_RTOL, err_msg=case["name"] + " ядро |S1|")
        np.testing.assert_allclose(r["core2"], r["ref2"], rtol=_RTOL, err_msg=case["name"] + " ядро |S2|")


def test_efficiencies_match_wiscombe():
    """Q_ext, Q_sca совпадают с таблицей Уискомба."""
    for case in _CASES:
        r = _evaluate(case)
        assert math.isclose(r["q_ext"], case["q_ext"], rel_tol=1e-5), case["name"] + " Qext"
        assert math.isclose(r["q_sca"], case["q_sca"], rel_tol=1e-5), case["name"] + " Qsca"


def test_absolute_rcs_matches_wiscombe():
    """Абсолютная σ(θ)=(4π/x²)|S|² ядра совпадает с σ из таблицы (R=1 м)."""
    for case in _CASES:
        r = _evaluate(case)
        th = np.deg2rad(_ANGLES.astype(float))
        sig = bistatic_rcs_sphere(1.0, [case["m"] ** 2], case["x"], th)
        pref = 4 * math.pi / case["x"] ** 2
        np.testing.assert_allclose(sig["sigma_E"], pref * r["ref2"] ** 2, rtol=_RTOL, err_msg=case["name"] + " σ_E")
        np.testing.assert_allclose(sig["sigma_H"], pref * r["ref1"] ** 2, rtol=_RTOL, err_msg=case["name"] + " σ_H")


if __name__ == "__main__":
    print("Верификация против эталона MIEV0 (Wiscombe 1980, Appl. Opt. 19:1505)\n")
    ok = True
    for case in _CASES:
        r = _evaluate(case)
        d1 = float(np.max(np.abs(r["core1"] - r["ref1"]) / r["ref1"]))
        d2 = float(np.max(np.abs(r["core2"] - r["ref2"]) / (r["ref2"] + 1e-30)))
        dq = max(abs(r["q_ext"] - case["q_ext"]) / case["q_ext"],
                 abs(r["q_sca"] - case["q_sca"]) / case["q_sca"])
        good = d1 < _RTOL and d2 < _RTOL and dq < 1e-5
        ok &= good
        print(f"  [{'OK' if good else 'FAIL'}] {case['name']:34s} "
              f"max|Δ|S1|={d1:.1e}  |S2|={d2:.1e}  Qext/Qsca={dq:.1e}")
    print("\nИТОГ:", "совпадение с опубликованной таблицей" if ok else "ЕСТЬ РАСХОЖДЕНИЯ")
    sys.exit(0 if ok else 1)
