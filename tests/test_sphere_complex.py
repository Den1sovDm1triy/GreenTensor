# SPDX-License-Identifier: MIT

"""Фаза 1 — проверки комплексной eps, поглощения и металл-диэлектрик (01_sphere.py).

Покрывает то, что раньше было сломано в обработке комплексной проницаемости:
  1. поглощающий диэлектрик: Q_sca и Q_abs совпадают с аналитическим Ми;
  2. знак поглощения: Im(eps)>0 -> Q_abs>0, Im(eps)<0 -> Q_abs<0 (конвенция e^{-iωt});
  3. реальный металл (Re(eps)<0): совпадение с точной mie_ab(sqrt(eps));
  4. PEC-предел (|eps| огромна): конечный результат, сходимость к mie_pec;
  5. переход металл-диэлектрик: металлическая оболочка экранирует ядро -> PEC сферы.

Запуск:
    python3 tests/test_sphere_complex.py
    pytest tests/test_sphere_complex.py
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from analytic_mie import mie_ab, mie_pec, _q_sca, _q_ext  # noqa: E402
from _loader import load_sphere  # noqa: E402

warnings.filterwarnings("ignore")
_sphere = load_sphere()
RCSCalculator = _sphere.RCSCalculator


def _coeffs(k0, a, eps, toch):
    n = len(a)
    c = RCSCalculator(k0=k0, toch=toch, n=n, phi=1.0, a=list(a), eps=list(eps),
                      miy=[1.0] * n, k1=k0)
    c.run_calculation()
    return np.asarray(c.Mn, complex), np.asarray(c.Nn, complex)


def _q(Mn, Nn, x):
    nn = np.arange(1, len(Mn) + 1)
    qs = float((2 / x**2) * np.sum((2 * nn + 1) * (np.abs(Mn) ** 2 + np.abs(Nn) ** 2)))
    qe = float((2 / x**2) * np.sum((2 * nn + 1) * np.real(Mn + Nn)))
    return qs, qe - qs  # (Q_sca, Q_abs)


def _homog(eps, x, toch):
    return _coeffs(x, [0.5, 1.0, 1.0], [eps, eps, 1.0], toch)


def test_absorption_matches_mie():
    """Поглощающий диэлектрик: Q_sca и Q_abs совпадают с аналитикой."""
    print("\n[1] Поглощение eps=2.25+i*loss (Q_sca и Q_abs vs Ми):")
    x, toch = 3.0, 12
    for loss in (0.01, 0.1, 0.5, 1.0, 5.0):
        eps = 2.25 + 1j * loss
        gsca, gabs = _q(*_homog(eps, x, toch), x)
        n, a, b = mie_ab(np.sqrt(eps), x, toch)
        rsca, rabs = _q_sca(n, a, b, x), _q_ext(n, a, b, x) - _q_sca(n, a, b, x)
        ds, da = abs(gsca - rsca) / rsca, abs(gabs - rabs) / abs(rabs)
        print(f"    loss={loss:<5} dQsca={ds:.1e} dQabs={da:.1e}")
        assert ds < 1e-6 and da < 1e-6, f"loss={loss}: Qsca/Qabs не совпали ({ds:.1e}/{da:.1e})"


def test_absorption_sign():
    """Im(eps)>0 поглощает (Q_abs>0); Im(eps)<0 — усиление (Q_abs<0)."""
    print("\n[2] Знак поглощения по знаку Im(eps):")
    x, toch = 3.0, 12
    _, qabs_pos = _q(*_homog(2.25 + 0.2j, x, toch), x)
    _, qabs_neg = _q(*_homog(2.25 - 0.2j, x, toch), x)
    print(f"    Im>0: Qabs={qabs_pos:+.3e}   Im<0: Qabs={qabs_neg:+.3e}")
    assert qabs_pos > 0, "Im(eps)>0 должно поглощать (Q_abs>0)"
    assert qabs_neg < 0, "Im(eps)<0 должно усиливать (Q_abs<0)"


def test_real_metal_matches_mie():
    """Реальный металл (Re(eps)<0): совпадение с точной mie_ab(sqrt(eps))."""
    print("\n[3] Реальный металл vs точная аналитика:")
    x, toch = 3.0, 12
    for eps in (-10 + 1j, -100 + 1j, -100 - 1j):
        gsca, _ = _q(*_homog(eps, x, toch), x)
        n, a, b = mie_ab(np.sqrt(eps), x, toch)
        rel = abs(gsca - _q_sca(n, a, b, x)) / _q_sca(n, a, b, x)
        print(f"    eps={eps!s:>10}: dQsca={rel:.1e}")
        assert rel < 1e-3, f"eps={eps}: металл не совпал с аналитикой ({rel:.1e})"


def test_pec_limit_is_finite_and_correct():
    """Огромная |eps| (PEC-прокси) даёт конечный результат и сходится к PEC."""
    print("\n[4] PEC-предел (|eps| огромна):")
    x, toch = 3.0, 12
    pec = _q_sca(*mie_pec(x, toch), x)
    for eps in (-1.7e7j, 1e7j, 1e8 + 1e8j):
        gsca, _ = _q(*_homog(eps, x, toch), x)
        assert np.isfinite(gsca), f"eps={eps}: переполнение (nan/inf)"
        rel = abs(gsca - pec) / pec
        print(f"    eps={eps!s:>14}: GT={gsca:.5e} PEC={pec:.5e} откл={rel:.1%}")
        assert rel < 2e-2, f"eps={eps}: не сошлось к PEC ({rel:.1%})"


def test_metal_shell_shields_to_pec():
    """Диэлектрическое ядро + металлическая оболочка => экранирование => PEC сферы r=1."""
    print("\n[5] Металлическая оболочка экранирует ядро (-> PEC):")
    x, toch = 3.0, 12
    pec = _q_sca(*mie_pec(x, toch), x)
    gsca, _ = _q(*_coeffs(x, [0.5, 1.0, 1.0], [4.0, -1.7e7j, 1.0], toch), x)
    rel = abs(gsca - pec) / pec
    print(f"    диэл-ядро+металл-оболочка: GT={gsca:.5e} PEC={pec:.5e} откл={rel:.1%}")
    assert rel < 2e-2, f"экранирование металлической оболочкой нарушено ({rel:.1%})"


_TESTS = [
    test_absorption_matches_mie,
    test_absorption_sign,
    test_real_metal_matches_mie,
    test_pec_limit_is_finite_and_correct,
    test_metal_shell_shields_to_pec,
]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  FAIL {fn.__name__}: {exc}")
            ok = False
    print("\nOK: Комплексные/металл проверки пройдены." if ok else "\nFAIL: Есть провалы.")
    sys.exit(0 if ok else 1)
