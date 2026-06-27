# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Magnetodielectric checks for the canonical 01_sphere.py solver.

The arbiter is the closed homogeneous-sphere Mie formula with independent eps
and mu. This file intentionally targets ``RCSCalculator`` from ``01_sphere.py``:
that module is the canonical layered-sphere solver in this checkout.
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _loader import load_sphere  # noqa: E402
from analytic_mie import _q_ext, _q_sca, mie_ab, mie_ab_eps_mu, mie_pec  # noqa: E402

warnings.filterwarnings("ignore")
_sphere = load_sphere()
RCSCalculator = _sphere.RCSCalculator


def _coeffs_layers(a, eps, miy, x: float, toch: int):
    calc = RCSCalculator(
        k0=x,
        toch=toch,
        n=len(a),
        phi=1.0,
        a=list(a),
        eps=list(eps),
        miy=list(miy),
        k1=x,
    )
    calc.run_calculation()
    return np.asarray(calc.Mn, complex), np.asarray(calc.Nn, complex)


def _coeffs(eps: complex, mu: complex, x: float, toch: int):
    return _coeffs_layers([0.5, 1.0, 1.0], [eps, eps, 1.0], [mu, mu, 1.0], x, toch)


def _q(Mn, Nn, x: float):
    n = np.arange(1, len(Mn) + 1)
    return _q_sca(n, Mn, Nn, x), _q_ext(n, Mn, Nn, x)


def test_eps_mu_formula_reduces_to_nonmagnetic_mie():
    n, a, b = mie_ab_eps_mu(2.25 + 0.1j, 1.0, 1.4, 8)
    n0, a0, b0 = mie_ab(np.sqrt(2.25 + 0.1j), 1.4, 8)
    assert np.array_equal(n, n0)
    assert np.allclose(a, a0, rtol=1e-12, atol=1e-12)
    assert np.allclose(b, b0, rtol=1e-12, atol=1e-12)


def test_homogeneous_magnetodielectric_coefficients_match_closed_mie():
    cases = [
        (1.0, 4.0, 0.5),
        (2.25, 1.7, 1.2),
        (2.25 + 0.08j, 1.4 + 0.03j, 1.0),
    ]
    for eps, mu, x in cases:
        toch = 10
        Mn, Nn = _coeffs(eps, mu, x, toch)
        n, a, b = mie_ab_eps_mu(eps, mu, x, toch)
        assert np.allclose(Mn, np.conj(a), rtol=2e-6, atol=2e-8), (eps, mu, x, Mn[:3], np.conj(a[:3]))
        assert np.allclose(Nn, np.conj(b), rtol=2e-6, atol=2e-8), (eps, mu, x, Nn[:3], np.conj(b[:3]))


def test_homogeneous_magnetodielectric_cross_sections_match_closed_mie():
    eps, mu, x, toch = 2.25 + 0.05j, 1.6 + 0.02j, 1.5, 12
    Mn, Nn = _coeffs(eps, mu, x, toch)
    n, a, b = mie_ab_eps_mu(eps, mu, x, toch)
    got_sca, got_ext = _q(Mn, Nn, x)
    ref_sca, ref_ext = _q_sca(n, a, b, x), _q_ext(n, a, b, x)
    assert abs(got_sca - ref_sca) / ref_sca < 2e-6
    assert abs(got_ext - ref_ext) / abs(ref_ext) < 2e-6


def test_mu_only_sphere_has_distinct_electric_and_magnetic_response():
    eps, mu, x, toch = 1.0, 4.0, 0.5, 10
    Mn, Nn = _coeffs(eps, mu, x, toch)
    _, a, b = mie_ab_eps_mu(eps, mu, x, toch)
    assert abs(np.conj(a[0]) - np.conj(b[0])) > 1e-2
    assert abs(Mn[0] - Nn[0]) > 1e-2


def test_dual_sphere_zero_backscatter_and_equal_coefficients():
    """For eps=mu the Kerker duality condition gives a_n=b_n and Q_back=0."""
    eps, mu, x, toch = 4.0, 4.0, 2.0, 14
    Mn, Nn = _coeffs(eps, mu, x, toch)
    n = np.arange(1, toch + 1)
    q_back = float((1.0 / x**2) * abs(np.sum((2 * n + 1) * ((-1) ** n) * (Mn - Nn))) ** 2)
    got_sca, got_ext = _q(Mn, Nn, x)

    assert np.allclose(Mn, Nn, rtol=1e-12, atol=1e-12)
    assert q_back < 1e-24
    assert abs(got_ext - got_sca) < 1e-12 * (1.0 + got_sca)


def test_lossless_magnetodielectric_energy_conservation():
    """For real positive eps,mu the passive lossless sphere has Q_ext=Q_sca."""
    eps, mu, x, toch = 2.25, 1.7, 1.2, 12
    Mn, Nn = _coeffs(eps, mu, x, toch)
    got_sca, got_ext = _q(Mn, Nn, x)

    assert got_sca > 0.0
    assert abs(got_ext - got_sca) < 1e-12 * (1.0 + got_sca)


def test_magnetodielectric_subdivision_invariance():
    eps, mu, x, toch = 2.4 + 0.04j, 1.6 + 0.02j, 1.3, 12
    one = _coeffs(eps, mu, x, toch)
    split = _coeffs_layers(
        [0.25, 0.50, 0.75, 1.0, 1.0],
        [eps, eps, eps, eps, 1.0],
        [mu, mu, mu, mu, 1.0],
        x,
        toch,
    )
    assert np.allclose(split[0], one[0], rtol=2e-6, atol=2e-8)
    assert np.allclose(split[1], one[1], rtol=2e-6, atol=2e-8)


def test_magnetic_core_metal_shell_shields_to_pec():
    x, toch = 3.0, 12
    Mn, Nn = _coeffs_layers(
        [0.5, 1.0, 1.0],
        [1.0, -1.7e7j, 1.0],
        [4.0, 1.0, 1.0],
        x,
        toch,
    )
    got_sca, _ = _q(Mn, Nn, x)
    ref_sca = _q_sca(*mie_pec(x, toch), x)
    assert abs(got_sca - ref_sca) / ref_sca < 2e-2


_TESTS = [
    test_eps_mu_formula_reduces_to_nonmagnetic_mie,
    test_homogeneous_magnetodielectric_coefficients_match_closed_mie,
    test_homogeneous_magnetodielectric_cross_sections_match_closed_mie,
    test_mu_only_sphere_has_distinct_electric_and_magnetic_response,
    test_dual_sphere_zero_backscatter_and_equal_coefficients,
    test_lossless_magnetodielectric_energy_conservation,
    test_magnetodielectric_subdivision_invariance,
    test_magnetic_core_metal_shell_shields_to_pec,
]


if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"FAIL {fn.__name__}: {exc}")
            ok = False
    print("OK: magnetic sphere checks passed." if ok else "FAIL: magnetic sphere checks failed.")
    sys.exit(0 if ok else 1)
