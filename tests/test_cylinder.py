# SPDX-License-Identifier: MIT

"""Проверки бесконечного цилиндра (green_tensor/cylinder.py).

  [1] сохранение энергии без потерь: Q_ext = Q_sca (TM и TE);
  [2] знак поглощения: Im(m)>0 ⇒ Q_abs>0, Q_ext>Q_sca;
  [3] положительность и убывание Q_sca при x→0 (малый цилиндр рассеивает слабо);
  [4] конечный цилиндр → честный NotImplementedError.

Запуск: python3 tests/test_cylinder.py | pytest tests/test_cylinder.py
"""
from __future__ import annotations

import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from green_tensor import cylinder as cyl  # noqa: E402


def test_energy_conservation_lossless():
    print("\n[1] сохранение энергии (lossless): Q_ext=Q_sca:")
    m = math.sqrt(2.25)
    for mode in ("TM", "TE"):
        for x in (0.5, 2.0, 5.0):
            cs = cyl.cross_sections_infinite(m, x, mode)
            d = abs(cs["q_ext"] - cs["q_sca"]) / cs["q_sca"]
            print(f"    {mode} x={x}: Qext={cs['q_ext']:.4f} Qsca={cs['q_sca']:.4f} |Qabs|/Qsca={d:.1e}")
            assert d < 1e-9, f"нарушение энергобаланса {mode} x={x}: {d:.1e}"


def test_absorption_sign():
    print("\n[2] поглощение Im(m)>0 ⇒ Q_abs>0:")
    m = complex(math.sqrt(2.25), 0.05)
    for mode in ("TM", "TE"):
        cs = cyl.cross_sections_infinite(m, 2.0, mode)
        print(f"    {mode}: Q_abs={cs['q_abs']:.4e}")
        assert cs["q_abs"] > 0 and cs["q_ext"] > cs["q_sca"]


def test_small_x_positive_decreasing():
    print("\n[3] малый цилиндр: Q_sca>0 и убывает с x→0:")
    m = math.sqrt(2.25)
    q = [cyl.cross_sections_infinite(m, x, "TM")["q_sca"] for x in (0.05, 0.1, 0.2)]
    print(f"    Q_sca(x=0.05,0.1,0.2) = {[f'{v:.2e}' for v in q]}")
    assert all(v > 0 for v in q) and q[0] < q[1] < q[2]


def test_finite_not_implemented():
    print("\n[4] конечный цилиндр → NotImplementedError:")
    try:
        cyl.finite()
        raised = False
    except NotImplementedError:
        raised = True
    assert raised
    print("    честный NotImplementedError — ок")


_TESTS = [test_energy_conservation_lossless, test_absorption_sign,
          test_small_x_positive_decreasing, test_finite_not_implemented]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ cylinder проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
