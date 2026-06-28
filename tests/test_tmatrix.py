# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Проверки адаптера T-матрицы (green_tensor/tmatrix.py).

Спина архитектуры: примитив -> T-матрица в сферическом базисе ВСВФ -> сечения.
Валидируется против НЕЗАВИСИМОГО аналитического арбитра Ми (Bohren & Huffman,
tests/analytic_mie.py), а также проверяется согласованность T-матрицы, собранной
из канонического ядра сферы (01_sphere/sphere_core), с прямой сборкой из коэффициентов a_n, b_n.

Запуск:
    python3 tests/test_tmatrix.py
    pytest tests/test_tmatrix.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from analytic_mie import mie_ab, q_sca, q_ext, wiscombe_nmax  # noqa: E402
from green_tensor import sphere_core as mc  # noqa: E402
from green_tensor import tmatrix as tm  # noqa: E402

# (eps слоя, толеранс) — диэлектрик, поглощающий, «металл» (PEC-прокси)
_CASES = [(2.25, 1e-9), (2.25 + 0.3j, 1e-9), (4.0 + 0.5j, 1e-9)]
_X = (0.5, 2.0, 6.0)


def _q_back_analytic(m, x, nmax):
    n, a, b = mie_ab(m, x, nmax)
    return float((1.0 / x**2) * np.abs(np.sum((2 * n + 1) * ((-1) ** n) * (a - b))) ** 2)


def test_sphere_tmatrix_cross_sections_vs_analytic():
    """Сечения из T-матрицы (адаптер sphere_core) == независимый Ми (B&H)."""
    print("\n[1] T-матрица сферы (адаптер sphere_core) vs аналитика B&H:")
    for eps, tol in _CASES:
        m = np.sqrt(eps)
        for x in _X:
            nmax = max(wiscombe_nmax(x), 6)
            sphere = mc.MieSphere(k0=x, a=[1.0], eps=[eps], toch=nmax)
            T = tm.sphere_tmatrix(sphere)
            d_sca = abs(T.q_sca() - q_sca(m, x, nmax)) / q_sca(m, x, nmax)
            d_ext = abs(T.q_ext() - q_ext(m, x, nmax)) / abs(q_ext(m, x, nmax))
            qb = _q_back_analytic(m, x, nmax)
            d_bk = abs(T.q_back() - qb) / qb
            print(f"    eps={eps!s:>12} x={x}: dQsca={d_sca:.1e} dQext={d_ext:.1e} dQback={d_bk:.1e}")
            assert d_sca < tol and d_ext < tol and d_bk < tol


def test_from_ab_optical_theorem():
    """Прямая сборка из a_n,b_n: оптическая теорема (Q_ext>=Q_sca, Q_abs знак)."""
    print("\n[2] from_ab + оптическая теорема:")
    x, nmax = 4.0, 12
    # без потерь: Q_abs ~ 0
    n, a, b = mie_ab(1.5, x, nmax)
    T = tm.from_ab(n, a, b, x)
    assert abs(T.q_abs()) < 1e-9, f"lossless Q_abs должно быть ~0, получено {T.q_abs():.1e}"
    assert abs(T.q_sca() - q_sca(1.5, x, nmax)) / q_sca(1.5, x, nmax) < 1e-12
    # поглощение: Q_abs > 0 и Q_ext > Q_sca
    n, a, b = mie_ab(1.5 + 0.05j, x, nmax)
    T = tm.from_ab(n, a, b, x)
    print(f"    Q_sca={T.q_sca():.4f} Q_ext={T.q_ext():.4f} Q_abs={T.q_abs():.4f}")
    assert T.q_abs() > 0 and T.q_ext() > T.q_sca()


def test_adapter_matches_direct_construction():
    """T(sphere_core) и T(analytic a,b) дают совпадающие сечения (согласованность конвенций)."""
    print("\n[3] согласованность sphere_core-адаптера и прямой сборки:")
    x, nmax = 3.0, 12
    for eps, _ in _CASES:
        T1 = tm.sphere_tmatrix(mc.MieSphere(k0=x, a=[1.0], eps=[eps], toch=nmax))
        n, a, b = mie_ab(np.sqrt(eps), x, nmax)
        T2 = tm.from_ab(n, a, b, x)
        d = abs(T1.q_sca() - T2.q_sca()) / T2.q_sca()
        print(f"    eps={eps!s:>12}: dQsca={d:.1e}")
        assert d < 1e-9


def test_tmatrix_validation():
    """Конструктор отвергает несогласованные/некорректные входы."""
    print("\n[4] валидация конструктора:")
    raised = 0
    bad = [
        dict(n=[1, 2], t_M=[0j], t_N=[0j, 0j], x=1.0),       # разные длины
        dict(n=[1], t_M=[0j], t_N=[0j], x=0.0),               # x не положителен
        dict(n=[1], t_M=[0j], t_N=[0j], x=float("inf")),      # x не конечен
    ]
    for cfg in bad:
        try:
            tm.DiagonalTMatrix(**cfg)
        except ValueError:
            raised += 1
    print(f"    отклонено: {raised}/3")
    assert raised == 3


_TESTS = [test_sphere_tmatrix_cross_sections_vs_analytic, test_from_ab_optical_theorem,
          test_adapter_matches_direct_construction, test_tmatrix_validation]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ tmatrix проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
