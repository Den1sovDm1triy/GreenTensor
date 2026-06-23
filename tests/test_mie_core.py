"""Фаза 2 — проверки чистого ядра green_tensor/mie_core.py.

Гарантируют, что вынесенное ядро:
  * совпадает с аналитическим Ми по сечениям (диэлектрик / поглощение / PEC);
  * совпадает с исправленным RCSCalculator по коэффициентам и линейной ДН;
  * поддерживает параметры polarization и problem;
  * валидирует входные параметры.

Запуск:
    python3 tests/test_mie_core.py
    pytest tests/test_mie_core.py
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from analytic_mie import (  # noqa: E402
    q_sca, mie_pec, _q_sca, mie_helicity, wiscombe_nmax, trapezoid,
)
from _loader import load_sphere  # noqa: E402
from green_tensor import mie_core as mc  # noqa: E402

warnings.filterwarnings("ignore")
_S = load_sphere()


def _rcs(k0, a, eps, toch, phi=1.0):
    c = _S.RCSCalculator(k0=k0, toch=toch, n=len(a), phi=phi, a=list(a),
                         eps=list(eps), miy=[1.0] * len(a), k1=k0)
    c.run_calculation()
    return np.asarray(c.Mn, complex), np.asarray(c.Nn, complex), c


def test_cross_sections_vs_mie():
    print("\n[1] mie_core сечения vs аналитика:")
    x, toch = 3.0, 12
    # диэлектрик и поглощающий — против точного Ми
    for eps in (2.25, 2.25 + 0.3j):
        cs = mc.MieSphere(k0=x, a=[0.5, 1, 1], eps=[eps, eps, 1.0], toch=toch).cross_sections()
        ref = q_sca(np.sqrt(eps), x, toch)
        rel = abs(cs["q_sca"] - ref) / ref
        print(f"    eps={eps!s:>14}: dQsca={rel:.1e}")
        assert rel < 1e-9, f"eps={eps}: Qsca не совпал ({rel:.1e})"
    # PEC-предел
    cs = mc.MieSphere(k0=x, a=[0.5, 1, 1], eps=[-1.7e7j, -1.7e7j, 1.0], toch=toch).cross_sections()
    pec = _q_sca(*mie_pec(x, toch), x)
    rel = abs(cs["q_sca"] - pec) / pec
    print(f"    PEC: dQsca={rel:.1e}")
    assert np.isfinite(cs["q_sca"]) and rel < 2e-2


def test_coefficients_match_rcs():
    print("\n[2] mie_core коэффициенты vs RCSCalculator:")
    x, toch = 3.0, 12
    for eps, tol in (([2.25, 2.25, 1.0], 1e-9),
                     ([2.25 + 0.3j, 2.25 + 0.3j, 1.0], 1e-9),
                     ([3.9, 2.5, 1.0], 1e-9),
                     ([-1.7e7j, -1.7e7j, 1.0], 5e-3)):
        Mn, Nn = mc.MieSphere(k0=x, a=[0.5, 1, 1], eps=eps, toch=toch).coefficients()
        rMn, rNn, _ = _rcs(x, [0.5, 1, 1], eps, toch)
        d = max(np.max(np.abs(Mn - rMn) / (np.abs(rMn) + 1e-30)),
                np.max(np.abs(Nn - rNn) / (np.abs(rNn) + 1e-30)))
        print(f"    eps={eps!s:>30}: max rel={d:.1e} (tol={tol:.0e})")
        assert d < tol, f"eps={eps}: коэффициенты разошлись с RCSCalculator ({d:.1e})"


def test_linear_matches_rcs():
    """Линейная ДН (оригинальная формула E_θ/E_φ) совпадает с RCSCalculator."""
    print("\n[3] линейная E_θ/E_φ vs RCSCalculator (точный порт формулы):")
    x, toch = 6.0, 16
    for phi in (0.5, 1.0):
        p = mc.MieSphere(k0=x, a=[0.5, 1, 1], eps=[3.9, 2.5, 1.0], toch=toch,
                         polarization="linear", phi=phi).pattern()
        _, _, c = _rcs(x, [0.5, 1, 1], [3.9, 2.5, 1.0], toch, phi=phi)
        Et, Ep = np.asarray(c.E_teta[0]), np.asarray(c.E_phi[0])
        d_t = np.max(np.abs(p["E_theta"] - Et)) / np.max(np.abs(Et))
        d_p = np.max(np.abs(p["E_phi"] - Ep)) / np.max(np.abs(Ep))
        print(f"    phi={phi}: d(E_θ)={d_t:.1e} d(E_φ)={d_p:.1e}")
        assert d_t < 1e-9 and d_p < 1e-9, f"phi={phi}: линейная ДН разошлась с RCSCalculator"


def test_polarization_and_problem_params():
    print("\n[4] параметры polarization / problem:")
    x, toch = 3.0, 12
    base = dict(k0=x, a=[0.5, 1, 1], eps=[2.25, 2.25, 1.0], toch=toch)
    lin_h = mc.MieSphere(**base, polarization="linear", orientation="horizontal").pattern()
    lin_v = mc.MieSphere(**base, polarization="linear", orientation="vertical").pattern()
    cir = mc.MieSphere(**base, polarization="circular").pattern()
    assert {"E_theta", "E_phi", "E"} <= set(lin_h) and {"E_op", "E_kp"} <= set(cir)
    assert np.all(np.isfinite(lin_h["E"])) and np.all(np.isfinite(cir["E_op"]))
    # orientation выбирает компоненту: horizontal -> E_phi, vertical -> E_theta
    assert np.allclose(lin_h["E"], lin_h["E_phi"]) and np.allclose(lin_v["E"], lin_v["E_theta"])
    diff_hv = np.max(np.abs(lin_h["E"] - lin_v["E"]))
    print(f"    горизонтальная(E_φ) vs вертикальная(E_θ): max|разн|={diff_hv:.3e}")
    assert diff_hv > 1e-9, "горизонтальная и вертикальная должны различаться"
    Mn_d, _ = mc.MieSphere(**base, problem="diffraction").coefficients()
    Mn_a, _ = mc.MieSphere(**base, problem="antenna").coefficients()
    diff = np.max(np.abs(Mn_d - Mn_a))
    print(f"    diffraction vs antenna: max|dMn|={diff:.3e} (должно быть >0)")
    assert diff > 1e-6, "antenna и diffraction дали одинаковые коэффициенты"
    assert np.all(np.isfinite(Mn_a))


def test_invalid_params():
    print("\n[5] валидация параметров:")
    cases = [
        dict(k0=3.0, a=[0.5, 1, 1], eps=[2.25, 2.25, 1.0], toch=8, polarization="bogus"),
        dict(k0=3.0, a=[0.5, 1, 1], eps=[2.25, 2.25, 1.0], toch=8, problem="bogus"),
        dict(k0=3.0, a=[0.5, 1, 1], eps=[2.25, 2.25, 1.0], toch=8, orientation="bogus"),
        dict(k0=3.0, a=[0.5, 1, 1], eps=[1.0], toch=8),  # len(eps) != len(a)
    ]
    raised = 0
    for cfg in cases:
        try:
            mc.MieSphere(**cfg)
        except ValueError:
            raised += 1
    print(f"    отклонено некорректных конфигураций: {raised}/4")
    assert raised == 4


def test_circular_vs_analytic_helicity():
    """Круговая (формула репо, θ от обратного направления) = 2·|S_co/S_cross|(π−θ) аналитики."""
    print("\n[6] круговые E_op/E_kp vs аналитика спиральности (угол θ -> π−θ):")
    th = np.linspace(0.02, np.pi - 0.02, 80)
    for eps_r in (2.25, 2.25 + 0.3j):
        for x in (2.0, 6.0):
            toch = max(wiscombe_nmax(x), 8)
            p = mc.MieSphere(k0=x, a=[0.5, 1, 1], eps=[eps_r, eps_r, 1.0],
                             toch=toch, polarization="circular").pattern(th)
            Sco, Scr = mie_helicity(np.sqrt(eps_r), x, np.pi - th, toch)  # развёрнутый угол
            dco = np.max(np.abs(p["E_op"] - 2 * np.abs(Sco))) / np.max(2 * np.abs(Sco))
            dcr = np.max(np.abs(p["E_kp"] - 2 * np.abs(Scr))) / np.max(2 * np.abs(Scr))
            print(f"    eps={eps_r!s:>12} x={x}: d(E_op)={dco:.1e} d(E_kp)={dcr:.1e}")
            assert dco < 1e-9 and dcr < 1e-9


def test_circular_integral_and_physics():
    """Круговая (в т.ч. многослойная): интеграл -> Q_sca; θ=0(назад) co->0, θ=π(вперёд) cross->0."""
    print("\n[7] круговая: интеграл Q_sca + forward/backscatter:")
    thF = np.linspace(1e-4, np.pi, 4000)
    for eps in ([2.25, 2.25, 1.0], [3.9, 2.5, 1.0], [2.25 + 0.3j, 2.25 + 0.3j, 1.0]):
        for x in (2.0, 6.0):
            toch = max(wiscombe_nmax(x), 8)
            m = mc.MieSphere(k0=x, a=[0.5, 1, 1], eps=eps, toch=toch, polarization="circular")
            pI = m.pattern(thF)
            # E_op=2|S_co(π−θ)|, E_kp=2|S_cross(π−θ)| => (1/(2x²))∫(E_op²+E_kp²)sinθ = Q_sca
            qI = (1 / (2 * x**2)) * trapezoid((pI["E_op"] ** 2 + pI["E_kp"] ** 2) * np.sin(thF), thF)
            q = m.cross_sections()["q_sca"]
            rel = abs(qI - q) / q
            pb = m.pattern(np.array([1e-3]))           # θ→0 = обратное направление
            pf = m.pattern(np.array([np.pi - 1e-3]))   # θ→π = прямое направление
            co_back = pb["E_op"][0] / (pb["E_kp"][0] + 1e-30)    # со/кросс в backscatter -> 0
            cross_fwd = pf["E_kp"][0] / (pf["E_op"][0] + 1e-30)  # кросс/со в forward -> 0
            print(f"    eps={eps!s:>26} x={x}: dQsca={rel:.1e} co_back={co_back:.1e} cross_fwd={cross_fwd:.1e}")
            assert rel < 1e-4, "интеграл круговой != Q_sca"
            assert co_back < 1e-3 and cross_fwd < 1e-3, "нарушены forward/backscatter свойства"


_TESTS = [test_cross_sections_vs_mie, test_coefficients_match_rcs,
          test_linear_matches_rcs, test_circular_vs_analytic_helicity,
          test_circular_integral_and_physics, test_polarization_and_problem_params,
          test_invalid_params]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ mie_core проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
