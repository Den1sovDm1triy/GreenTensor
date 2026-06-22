"""Фаза 1 — арбитр корректности для рассеяния на сфере (01_sphere.py).

Сравнивает сечение рассеяния Q_sca, вычисленное RCSCalculator для ОДНОРОДНОЙ
диэлектрической сферы, с независимым аналитическим рядом Ми (analytic_mie.py).

Запуск:
    python3 tests/test_sphere_mie.py        # standalone, без pytest
    pytest tests/test_sphere_mie.py         # если установлен pytest
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analytic_mie import q_sca, wiscombe_nmax  # noqa: E402
from _loader import load_sphere  # noqa: E402

_sphere = load_sphere()
RCSCalculator = _sphere.RCSCalculator

# Однородная диэлектрическая сфера: eps_r = 2.25 (m = 1.5), радиус нормирован к 1.
EPS_R = 2.25
M = np.sqrt(EPS_R)
X_VALUES = (0.5, 1.0, 2.0, 3.0, 5.0)
RTOL = 2e-2  # 2 %


def gt_q_sca(x: float, toch: int) -> float:
    """Q_sca от RCSCalculator для однородной сферы.

    Однослойную сферу моделируем ДВУМЯ идентичными слоями (eps одинаков) + воздух
    — корректный многослойный решатель обязан давать тот же ответ, что и одиночный
    слой (инвариантность к разбиению однородной области).
    """
    calc = RCSCalculator(
        k0=x, toch=toch, n=3, phi=1.0,
        a=[0.5, 1.0, 1.0], eps=[EPS_R, EPS_R, 1.0], miy=[1.0, 1.0, 1.0], k1=x,
    )
    calc.run_calculation()
    Mn = np.asarray(calc.Mn, dtype=complex)
    Nn = np.asarray(calc.Nn, dtype=complex)
    nn = np.arange(1, len(Mn) + 1)
    return float((2.0 / x**2) * np.sum((2 * nn + 1) * (np.abs(Mn) ** 2 + np.abs(Nn) ** 2)))


def _compare():
    rows = []
    for x in X_VALUES:
        toch = max(wiscombe_nmax(x), 5)
        ref = q_sca(M, x)
        try:
            got = gt_q_sca(x, toch)
            rel = abs(got - ref) / ref
            rows.append((x, toch, ref, got, rel, None))
        except Exception as exc:  # noqa: BLE001
            rows.append((x, toch, ref, None, None, repr(exc)))
    return rows


def test_sphere_qsca_matches_analytic_mie():
    rows = _compare()
    print(f"\nОднородная сфера, eps_r={EPS_R} (m={M:.3f}); сравнение Q_sca:")
    print(f"{'x':>5} {'toch':>5} {'Mie(эталон)':>14} {'GreenTensor':>14} {'отн.ошибка':>11}")
    worst = 0.0
    errors = []
    for x, toch, ref, got, rel, err in rows:
        if err is not None:
            print(f"{x:>5} {toch:>5} {ref:>14.6e} {'ОШИБКА':>14} {'-':>11}  {err}")
            errors.append((x, err))
            continue
        flag = "  <== расхождение" if rel > RTOL else ""
        print(f"{x:>5} {toch:>5} {ref:>14.6e} {got:>14.6e} {rel:>10.2%}{flag}")
        worst = max(worst, rel)
    if errors:
        raise AssertionError(f"RCSCalculator упал на {len(errors)} конфигурациях: {errors}")
    assert worst <= RTOL, f"Максимальное расхождение {worst:.2%} > допуска {RTOL:.0%}"


if __name__ == "__main__":
    try:
        test_sphere_qsca_matches_analytic_mie()
        print("\n✅ PASS: GreenTensor совпадает с аналитическим Ми в пределах допуска.")
    except AssertionError as exc:
        print(f"\n❌ FAIL: {exc}")
        sys.exit(1)
