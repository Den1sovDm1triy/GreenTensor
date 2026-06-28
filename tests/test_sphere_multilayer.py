# SPDX-License-Identifier: MIT

"""Фаза 1 — проверки многослойной рекурсии 01_sphere.py против аналитики Ми.

Три конвенционно-устойчивых теста, не требующих внешних таблиц:
  1. инвариантность к разбиению однородной сферы на N идентичных слоёв;
  2. «вакуумная оболочка» (прозрачный внешний слой) = меньшая однородная сфера;
  3. обратное рассеяние Q_b (чувствительно к знаку/конвенции) = аналитический Ми.

Запуск:
    python3 tests/test_sphere_multilayer.py    # standalone
    pytest tests/test_sphere_multilayer.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analytic_mie import mie_ab, q_sca, wiscombe_nmax  # noqa: E402
from _loader import load_sphere  # noqa: E402

_sphere = load_sphere()
RCSCalculator = _sphere.RCSCalculator

EPS_R = 2.25
M = np.sqrt(EPS_R)
RTOL = 1e-6


def _run(k0, a, eps, toch):
    n = len(a)
    calc = RCSCalculator(k0=k0, toch=toch, n=n, phi=1.0, a=list(a), eps=list(eps),
                         miy=[1.0] * n, k1=k0)
    calc.run_calculation()
    return np.asarray(calc.Mn, complex), np.asarray(calc.Nn, complex)


def _coeff_sum(Mn, Nn):
    nn = np.arange(1, len(Mn) + 1)
    return float(np.sum((2 * nn + 1) * (np.abs(Mn) ** 2 + np.abs(Nn) ** 2)))


def test_subdivision_invariance():
    """Однородную сферу можно разбить на любое число идентичных слоёв."""
    print("\n[1] Инвариантность к разбиению (eps=2.25):")
    for x in (0.5, 1.0, 2.0, 3.0, 5.0):
        toch = max(wiscombe_nmax(x), 6)
        sums = []
        for n_shell in (2, 3, 4, 6):
            radii = [(i + 1) / n_shell for i in range(n_shell)]
            a = radii + [1.0]
            eps = [EPS_R] * n_shell + [1.0]
            sums.append(_coeff_sum(*_run(x, a, eps, toch)))
        spread = (max(sums) - min(sums)) / np.mean(sums)
        print(f"    x={x:<4} разброс по слоям = {spread:.2e}")
        assert spread < RTOL, f"x={x}: результат зависит от числа идентичных слоёв ({spread:.2e})"


def test_vacuum_shell_equals_core():
    """Прозрачная (вакуумная) внешняя оболочка не должна менять рассеяние ядра."""
    print("\n[2] Вакуумная оболочка = меньшая однородная сфера:")
    for k0 in (2.0, 4.0, 6.0, 10.0):
        toch = max(wiscombe_nmax(k0), 8)
        got = _coeff_sum(*_run(k0, a=[0.5, 1.0, 1.0], eps=[EPS_R, 1.0, 1.0], toch=toch))
        xp = k0 * 0.5
        n, a, b = mie_ab(M, xp, toch)
        ref = float(np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2)))
        rel = abs(got - ref) / ref
        print(f"    k0={k0:<4} (x_ядра={xp}): откл = {rel:.2e}")
        assert rel < RTOL, f"k0={k0}: вакуумная оболочка искажает результат ({rel:.2e})"


def test_backscatter_matches_mie():
    """Обратное рассеяние Q_b — чувствительно к знаку/конвенции."""
    print("\n[3] Обратное рассеяние Q_b vs аналитический Ми:")
    for x in (0.5, 1.0, 2.0, 3.0, 5.0):
        toch = max(wiscombe_nmax(x), 6)
        Mn, Nn = _run(x, a=[0.5, 1.0, 1.0], eps=[EPS_R, EPS_R, 1.0], toch=toch)
        nn = np.arange(1, len(Mn) + 1)
        got = float((1 / x**2) * np.abs(np.sum((2 * nn + 1) * ((-1) ** nn) * (Mn - Nn))) ** 2)
        n, a, b = mie_ab(M, x, toch)
        ref = float((1 / x**2) * np.abs(np.sum((2 * n + 1) * ((-1) ** n) * (a - b))) ** 2)
        rel = abs(got - ref) / ref
        print(f"    x={x:<4}: откл = {rel:.2e}")
        assert rel < 1e-4, f"x={x}: обратное рассеяние не совпало ({rel:.2e})"


if __name__ == "__main__":
    ok = True
    for fn in (test_subdivision_invariance, test_vacuum_shell_equals_core,
               test_backscatter_matches_mie):
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ Все многослойные проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
