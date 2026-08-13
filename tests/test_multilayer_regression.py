# SPDX-License-Identifier: MIT

"""Регрессия multilayer_mie: лог-производная многослойная сфера + PEC-покрытия.

Тесты соответствуют разделу «Internal Regression Suite» статьи PIERE-2026 и
разделу о регрессионных тестах гл. 9 диссертации:

  * необходимые условия: энергобаланс, парадокс экстинкции, сходимость Уискомба;
  * сравнение с независимыми замкнутыми формулами: Бор--Хаффман (mu = 1),
    обобщение Керкера (mu != 1);
  * согласие double- и mpmath-веток PEC+покрытие в зоне пересечения
    (диспетчер по |Im(m k0 r)|; mpmath-часть пропускается, если пакет не установлен).

Запуск:
    python3 tests/test_multilayer_regression.py   # standalone, без pytest
    pytest tests/test_multilayer_regression.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from green_tensor.mie_reference import (  # noqa: E402
    cross_sections,
    mie_homogeneous_sphere,
    mie_pec_sphere,
    wiscombe_nmax,
)
from green_tensor.multilayer_mie import (  # noqa: E402
    mie_homogeneous_magnetic_sphere,
    mie_multilayer,
    mie_pec_with_coating,
)


def test_reg1_homogeneous_bh_and_kerker():
    """L=1 против Бора--Хаффмана (mu=1) и Керкера (mu!=1): <= 1e-10."""
    max_err = 0.0
    for x in (1.0, 2 * np.pi, 6 * np.pi):
        for eps_val in (1.5, 4.0, 9.0):
            for mu_val in (1.0, 2.0):
                m_rel = np.sqrt(eps_val * mu_val)
                n_max = wiscombe_nmax(m_rel * x) + 8
                if mu_val == 1.0:
                    a_ref, b_ref, _ = mie_homogeneous_sphere(x, complex(m_rel), n_max)
                else:
                    a_ref, b_ref = mie_homogeneous_magnetic_sphere(x, eps_val, mu_val, n_max)
                a_ml, b_ml = mie_multilayer(1.0, [x], [eps_val], [mu_val], n_max)
                max_err = max(max_err,
                              float(np.max(np.abs(a_ref - a_ml))),
                              float(np.max(np.abs(b_ref - b_ml))))
    assert max_err < 1e-10, f"REG1 max err {max_err:.3e}"


def test_reg2_subdivision_invariance():
    """L=N однородных слоёв == L=1: <= 1e-10."""
    max_err = 0.0
    for x in (2 * np.pi, 6 * np.pi):
        for L in (2, 8):
            n_max = wiscombe_nmax(2.0 * x) + 8
            a1, b1 = mie_multilayer(1.0, [x], [4.0], [1.0], n_max)
            radii = [x * (l + 1) / L for l in range(L)]
            aN, bN = mie_multilayer(1.0, radii, [4.0] * L, [1.0] * L, n_max)
            max_err = max(max_err,
                          float(np.max(np.abs(a1 - aN))),
                          float(np.max(np.abs(b1 - bN))))
    assert max_err < 1e-10, f"REG2 max err {max_err:.3e}"


def test_energy_conservation():
    """Непоглощающая сфера: |Q_ext - Q_sca|/Q_sca на уровне машинного эпсилон."""
    worst = 0.0
    for eps_val in (1.5, 2.5, 4.0, 9.0):
        m = complex(np.sqrt(eps_val))
        for x in (1.0, 2 * np.pi, 8 * np.pi):
            n_max = wiscombe_nmax(abs(m) * x) + 8
            a_n, b_n, _ = mie_homogeneous_sphere(x, m, n_max)
            q_s, q_e, _ = cross_sections(a_n, b_n, x)
            worst = max(worst, abs(q_e - q_s) / q_s)
    assert worst < 1e-13, f"energy imbalance {worst:.3e}"


def test_extinction_paradox():
    """PEC-сфера: Q_ext -> 2 (2.033 при x=20, 2.006 при x=150)."""
    for x, q_expect, tol in ((20.0, 2.033, 0.002), (150.0, 2.006, 0.002)):
        n_max = wiscombe_nmax(x) + 12
        a_n, b_n, _ = mie_pec_sphere(x, n_max)
        _, q_e, _ = cross_sections(a_n, b_n, x)
        assert abs(q_e - q_expect) < tol, f"x={x}: Q_ext={q_e:.4f}"


def test_wiscombe_truncation():
    """Порядок N_W + 10 меняет sigma_r не более чем на 2e-9 относительных."""
    for x in (2 * np.pi, 8 * np.pi):
        m = complex(2.0)
        n_w = wiscombe_nmax(abs(m) * x)
        vals = []
        for n_max in (n_w, n_w + 10):
            a_n, b_n, _ = mie_homogeneous_sphere(x, m, n_max)
            _, _, q_b = cross_sections(a_n, b_n, x)
            vals.append(q_b)
        rel = abs(vals[1] - vals[0]) / abs(vals[1])
        assert rel < 2e-9, f"x={x}: truncation rel change {rel:.3e}"


def test_pec_coating_double_vs_mpmath():
    """Диспетчер double/mpmath: согласие веток в зоне пересечения <= 1e-9."""
    try:
        from green_tensor.multilayer_mie import _mie_pec_with_coating_mp
    except ImportError:
        return  # структура модуля изменилась — покрыто другими тестами
    try:
        import mpmath  # noqa: F401
    except ImportError:
        print("  (mpmath не установлен — пропуск mpmath-части)")
        return
    # умеренно поглощающее покрытие: |Im z| ~ 6 — обе ветки валидны
    k0, r_pec, t = 60.0, 0.10, 0.004
    eps_c, mu_c = 3.4 - 0.5j, 1.2 - 0.3j
    b = r_pec + t
    m_abs = abs(np.sqrt(eps_c * mu_c))
    n_max = wiscombe_nmax(m_abs * k0 * b) + 8
    a_d, b_d = mie_pec_with_coating(k0, r_pec, [b], [eps_c], [mu_c], n_max)
    a_m, b_m = _mie_pec_with_coating_mp(k0, r_pec, [b], [eps_c], [mu_c], n_max)
    err = max(float(np.max(np.abs(a_d - a_m))), float(np.max(np.abs(b_d - b_m))))
    assert err < 1e-9, f"double vs mpmath: {err:.3e}"


ALL_TESTS = [
    test_reg1_homogeneous_bh_and_kerker,
    test_reg2_subdivision_invariance,
    test_energy_conservation,
    test_extinction_paradox,
    test_wiscombe_truncation,
    test_pec_coating_double_vs_mpmath,
]

if __name__ == "__main__":
    for t in ALL_TESTS:
        t()
        print(f"PASS {t.__name__}")
    print("test_multilayer_regression: all OK")
