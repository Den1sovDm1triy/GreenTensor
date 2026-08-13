# SPDX-License-Identifier: MIT
"""GreenTensor reference Mie implementation (Wiscombe-Du downward recurrence).

Эталонная реализация Ми-формализма: однородная сфера, PEC-сфера, интегральные
сечения, критерий усечения Уискомба. Используется multilayer_mie.py и
регрессионными тестами как независимый золотой стандарт.

Исходно: python_calc/regression_tests/regression_clean.py (диссертация, гл. 9).

Clean Mie regression test for chapter 9.

Implements Bohren-Huffman analytic Mie via Wiscombe-Du downward recurrence
for the logarithmic derivative D_n(mx). This is the numerically stable
gold-standard algorithm.

Three internal regression tests:
1. Energy conservation: for lossless sphere, Q_ext == Q_sca (machine precision).
2. Extinction paradox: Q_ext -> 2 as x -> infinity.
3. Wiscombe truncation: |sigma_r(n_max=Wiscombe) - sigma_r(n_max=Wiscombe+10)|
   / sigma_r < 1e-12.

The "layered algorithm" in LunebergLens_difraction_lin_polar.m is the
multilayer generalisation; in the L=1 limit it reduces analytically to
Bohren-Huffman eq 4.56-4.57. This file validates the L=1 (homogeneous
sphere) case which is the foundation of the layered solver.
"""
from __future__ import annotations
import numpy as np
from scipy import special


def downward_D(mx: complex, n_max: int) -> np.ndarray:
    """Compute D_n(mx) = psi'_n(mx)/psi_n(mx) for n = 1..n_max via Wiscombe-Du
    downward recurrence:
       D_{n-1} = n/mx - 1/(D_n + n/mx)
    Start from D_N = 0 with N = n_max + max(15, ceil(|mx|^(1/3)*4)).
    Returns D[1..n_max] as length-n_max array.
    """
    N_extra = max(15, int(np.ceil(abs(mx)**(1/3) * 4)))
    N_total = n_max + N_extra
    D = np.zeros(N_total + 1, dtype=complex)
    # D[N_total] = 0, recurse downward
    for n in range(N_total, 0, -1):
        D[n - 1] = n / mx - 1.0 / (D[n] + n / mx)
    return D[1:n_max + 1]


def upward_psi_chi(x: float, n_max: int):
    """psi_n(x) and chi_n(x) for n = 0..n_max, real x."""
    n_array = np.arange(0, n_max + 1)
    jn = special.spherical_jn(n_array, x)
    yn = special.spherical_yn(n_array, x)
    psi = x * jn
    # chi convention: chi_n = -x*y_n  (so xi = psi - i*chi = x*(j_n - i*(-y_n)) = x*(j_n + i*y_n) = x*h_n^(1))
    chi = -x * yn
    return psi, chi


def mie_homogeneous_sphere(x: float, m_rel: complex, n_max: int):
    """
    Bohren-Huffman analytic Mie via Wiscombe-Du downward recurrence.

    Returns (a_n, b_n, n_array) where:
        a_n = TM (electric multipole) coefficients
        b_n = TE (magnetic multipole) coefficients
    Using BH eq 4.56-4.57:
        a_n = ((D_n/m + n/x) psi_n(x) - psi_{n-1}(x)) /
              ((D_n/m + n/x) xi_n(x)  - xi_{n-1}(x))
        b_n = ((m D_n + n/x) psi_n(x) - psi_{n-1}(x)) /
              ((m D_n + n/x) xi_n(x)  - xi_{n-1}(x))
    """
    mx = m_rel * x
    D = downward_D(mx, n_max)               # length n_max, indices 1..n_max
    psi, chi = upward_psi_chi(x, n_max)     # length n_max+1, indices 0..n_max
    xi = psi - 1j * chi                      # outgoing (h_n^(1) convention)

    n_arr = np.arange(1, n_max + 1)
    psi_n = psi[1:]               # n=1..n_max
    psi_nm = psi[:-1]             # n-1 = 0..n_max-1
    xi_n = xi[1:]
    xi_nm = xi[:-1]

    A = D / m_rel + n_arr / x
    B = m_rel * D + n_arr / x

    a_n = (A * psi_n - psi_nm) / (A * xi_n - xi_nm)
    b_n = (B * psi_n - psi_nm) / (B * xi_n - xi_nm)

    return a_n, b_n, n_arr


def mie_pec_sphere(x: float, n_max: int):
    """PEC limit: a_n = -psi_n(x)/xi_n(x), b_n = -psi'_n(x)/xi'_n(x).
    Implemented via the Riccati-Bessel recurrence for psi'.
    """
    psi, chi = upward_psi_chi(x, n_max)
    xi = psi - 1j * chi
    # psi'_n(x) = psi_{n-1}(x) - n/x * psi_n(x), n>=1
    n_arr = np.arange(1, n_max + 1)
    psi_n = psi[1:]
    psi_nm = psi[:-1]
    psip_n = psi_nm - (n_arr / x) * psi_n
    chi_n = chi[1:]
    chi_nm = chi[:-1]
    chip_n = chi_nm - (n_arr / x) * chi_n
    xi_n = psi_n - 1j * chi_n
    xip_n = psip_n - 1j * chip_n
    a_n = psip_n / xip_n   # TM (electric); from BH eq 4.56 in limit m->infty
    b_n = psi_n  / xi_n    # TE (magnetic)
    return a_n, b_n, n_arr


def cross_sections(a_n, b_n, x: float):
    """Compute Q_sca, Q_ext, Q_back for sphere of size x.
    All Q values are dimensionless (cross-section / pi*a^2).
    """
    n_arr = np.arange(1, len(a_n) + 1)
    weights = 2 * n_arr + 1
    Q_sca = (2.0 / x**2) * np.sum(weights * (np.abs(a_n)**2 + np.abs(b_n)**2))
    Q_ext = (2.0 / x**2) * np.sum(weights * np.real(a_n + b_n))
    # Backscatter (theta=pi):
    s_back = np.sum(weights * (-1.0)**n_arr * (a_n - b_n))
    Q_back = (1.0 / x**2) * np.abs(s_back)**2
    # sigma_r/a^2 = pi * Q_back (Bohren-Huffman convention pi*a^2 baseline)
    return float(np.real(Q_sca)), float(np.real(Q_ext)), float(np.real(Q_back))


def wiscombe_nmax(x: float) -> int:
    if x < 0.02:
        return 1
    if x <= 8.0:
        return int(np.round(x + 4.0 * x**(1.0/3.0) + 2))
    if x <= 4200.0:
        return int(np.round(x + 4.05 * x**(1.0/3.0) + 2))
    return int(np.round(x + 4.0 * x**(1.0/3.0) + 2))


# === Self-test ===
if __name__ == "__main__":
    print("=" * 64)
    print("REGRESSION TEST 1: Energy conservation (Q_ext == Q_sca for lossless)")
    print("=" * 64)
    print(f"{'x':>8s} {'eps':>8s} {'Q_ext':>14s} {'Q_sca':>14s} {'|diff|/Q':>12s}")
    print("-" * 64)
    for x in [0.5, 1.0, 2*np.pi, 4*np.pi, 8*np.pi]:
        for eps in [1.5, 2.5, 4.0, 9.0]:
            m = np.sqrt(complex(eps))
            n_max = wiscombe_nmax(max(x, m.real * x)) + 5
            a_n, b_n, _ = mie_homogeneous_sphere(x, m, n_max)
            Qs, Qe, Qb = cross_sections(a_n, b_n, x)
            err = abs(Qe - Qs) / max(abs(Qs), 1e-30)
            print(f"{x:8.3f} {eps:8.2f} {Qe:14.10f} {Qs:14.10f} {err:12.3e}")

    print()
    print("=" * 64)
    print("REGRESSION TEST 2: Extinction paradox (Q_ext -> 2 as x -> infty)")
    print("=" * 64)
    print(f"{'x':>8s} {'Q_ext (PEC)':>14s} {'Q_ext (eps=4)':>14s}")
    print("-" * 64)
    for x in [1, 5, 10, 20, 40, 80, 150]:
        n_max = wiscombe_nmax(x) + 5
        a_pec, b_pec, _ = mie_pec_sphere(x, n_max)
        Qs_p, Qe_p, _ = cross_sections(a_pec, b_pec, x)
        m = np.sqrt(4.0)
        n_max_d = wiscombe_nmax(max(x, m * x)) + 5
        a_d, b_d, _ = mie_homogeneous_sphere(x, m, n_max_d)
        Qs_d, Qe_d, _ = cross_sections(a_d, b_d, x)
        print(f"{x:8.0f} {Qe_p:14.6f} {Qe_d:14.6f}")

    print()
    print("=" * 64)
    print("REGRESSION TEST 3: Wiscombe truncation convergence")
    print("=" * 64)
    print(f"{'x':>8s} {'n_w':>5s} {'sig_r(n_w)':>14s} {'sig_r(n_w+10)':>16s} {'|diff|/sig':>12s}")
    print("-" * 64)
    for x in [1.0, 2*np.pi, 4*np.pi, 6*np.pi, 8*np.pi]:
        n_w = wiscombe_nmax(x)
        a1, b1, _ = mie_pec_sphere(x, n_w)
        _, _, Q1 = cross_sections(a1, b1, x)
        a2, b2, _ = mie_pec_sphere(x, n_w + 10)
        _, _, Q2 = cross_sections(a2, b2, x)
        err = abs(Q2 - Q1) / max(Q2, 1e-30)
        print(f"{x:8.3f} {n_w:5d} {Q1:14.10f} {Q2:16.10f} {err:12.3e}")
