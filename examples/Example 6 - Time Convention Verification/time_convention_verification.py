"""
Example 6 — Time-convention verification: e^{+j omega t} canon.

Verifies that the GreenTensor scattering coefficients M_n, N_n
(green_tensor/calc.py, RCSCalculator) are expressed in the dissertation
canon e^{+j omega t}:

    * outgoing radial dependence  ~ Riccati-Hankel of the SECOND kind
      zeta_n(z) = z * h_n^{(2)}(z) = psi_n(z) + j * chi_n(z),
    * lossy media carry eps = eps' - j eps''.

Implementation note. Internally calc.py evaluates the ratio with
scipy.special.hankel1 (h^{(1)}, native to the e^{-i omega t} convention)
and then applies an explicit conjugation

    Mn = Mn.real - Mn.imag * 1j        # calc.py, calculate_scattering_coefficients

which is mathematically identical to evaluating the same ratio with
h^{(2)} directly, because for real arguments conj(h^{(1)}) = h^{(2)} and
psi, chi, Z are real. This script PROVES that identity numerically:

  V1  single-layer sphere  (eps=2.5,  k0a = 2*pi):
        GreenTensor Mn, Nn  ==  a_n, b_n computed INDEPENDENTLY with
        zeta = psi + j*chi (direct h^{(2)}, no conjugation trick)
  V2  8-layer Fuchs-optimized Luneburg lens (k0a = 4*pi): same identity
      via a log-derivative multilayer solver written directly in the
      e^{+j omega t} convention
  V3  outgoing-wave check: d/dr arg[ zeta_n(k0 r) ] ~ -k0  < 0,
      i.e. together with e^{+j omega t} the phase fronts move OUTWARD;
      the h^{(1)} choice would give +k0 (an incoming wave in this canon)
  V4  energy conservation for the lossless lens:
        Re(Mn) = |Mn|^2  (unitarity circle |w - 1/2| = 1/2)
        sigma_ext = sigma_sca (optical theorem)
  V5  observable invariance: backscatter sweep sigma_r(k0a) of the
      8-layer lens, GreenTensor vs the canonical h^{(2)} solver

Run:
    python3 time_convention_verification.py

Reference: dissertation eq. (M_n, N_n from Z, Y) — Riccati-Hankel
h^{(2)} form; GreenTensor_Theory.pdf, sec. "Scattering coefficients".
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy import special

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- import GreenTensor from the repository root -------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, _ROOT)

from green_tensor.calc import RCSCalculator  # noqa: E402

TOL = 1e-12


# ==========================================================================
# Canonical e^{+j omega t} Mie solver — written DIRECTLY with h^{(2)}.
# No conjugation anywhere: the convention enters once, through
# zeta = psi + j*chi  (zeta_n(z) = z h_n^{(2)}(z), BH sign chi_n = -z y_n).
# ==========================================================================
def _psi(n_max: int, z: float) -> np.ndarray:
    n = np.arange(1, n_max + 1)
    return z * special.spherical_jn(n, z)


def _chi(n_max: int, z: float) -> np.ndarray:
    n = np.arange(1, n_max + 1)
    return -z * special.spherical_yn(n, z)


def _psip(n_max: int, z: float) -> np.ndarray:
    full = np.empty(n_max + 1)
    full[0] = math.sin(z)
    full[1:] = _psi(n_max, z)
    n = np.arange(1, n_max + 1)
    return full[n - 1] - (n / z) * full[n]


def _chip(n_max: int, z: float) -> np.ndarray:
    full = np.empty(n_max + 1)
    full[0] = math.cos(z)
    full[1:] = _chi(n_max, z)
    n = np.arange(1, n_max + 1)
    return full[n - 1] - (n / z) * full[n]


def _zeta_zetap(n_max: int, z: float):
    """zeta_n(z) = z h_n^{(2)}(z) = psi + j chi — the e^{+j omega t} outgoing wave."""
    return (_psi(n_max, z) + 1j * _chi(n_max, z),
            _psip(n_max, z) + 1j * _chip(n_max, z))


def mie_canon_multilayer(k0: float, radii, eps_arr, mu_arr, n_max: int):
    """Multilayer Mie a_n, b_n in the e^{+j omega t} canon (log-derivative).

    Lossless layers only (real eps, mu) — sufficient for the Luneburg-lens
    verification cases. The convention enters solely via zeta = psi + j chi
    in the final ratio; the internal D-propagation is real.
    """
    L = len(radii)
    eps = [float(e) for e in eps_arr]
    mu = [float(m) for m in mu_arr]
    m_l = [math.sqrt(eps[l] * mu[l]) for l in range(L)]

    z1 = m_l[0] * k0 * radii[0]
    D_TM = _psip(n_max, z1) / _psi(n_max, z1)
    D_TE = D_TM.copy()

    for l in range(1, L):
        z_in = m_l[l] * k0 * radii[l - 1]
        z_out = m_l[l] * k0 * radii[l]
        a_TM = (m_l[l - 1] * eps[l]) / (m_l[l] * eps[l - 1])
        a_TE = (m_l[l - 1] * mu[l]) / (m_l[l] * mu[l - 1])
        D_TM_in, D_TE_in = a_TM * D_TM, a_TE * D_TE

        psi_i, psip_i = _psi(n_max, z_in), _psip(n_max, z_in)
        chi_i, chip_i = _chi(n_max, z_in), _chip(n_max, z_in)
        psi_o, psip_o = _psi(n_max, z_out), _psip(n_max, z_out)
        chi_o, chip_o = _chi(n_max, z_out), _chip(n_max, z_out)
        alpha = psi_o * chi_i - psi_i * chi_o
        alpha_p = psip_o * chi_i - psi_i * chip_o
        gamma = psip_i * chi_o - psi_o * chip_i
        gamma_p = psip_i * chip_o - psip_o * chip_i
        D_TM = (gamma_p + alpha_p * D_TM_in) / (gamma + alpha * D_TM_in)
        D_TE = (gamma_p + alpha_p * D_TE_in) / (gamma + alpha * D_TE_in)

    D_TM_f = (m_l[-1] / eps[-1]) * D_TM
    D_TE_f = (m_l[-1] / mu[-1]) * D_TE

    y = k0 * radii[-1]
    psi_y, psip_y = _psi(n_max, y), _psip(n_max, y)
    zeta_y, zetap_y = _zeta_zetap(n_max, y)
    a_n = (psip_y - D_TM_f * psi_y) / (zetap_y - D_TM_f * zeta_y)
    b_n = (psip_y - D_TE_f * psi_y) / (zetap_y - D_TE_f * zeta_y)
    return a_n, b_n


# ==========================================================================
# Test configurations
# ==========================================================================
# 8-layer Fuchs-optimized Luneburg lens (dissertation, ch. 4 / B1 case)
FUCHS_RADII = [0.34, 0.49, 0.59, 0.69, 0.77, 0.84, 0.91, 0.97]
FUCHS_EPS = [1.94, 1.82, 1.71, 1.59, 1.47, 1.35, 1.24, 1.12]


def gt_coefficients(k0: float, radii, eps_arr, n_terms: int):
    """Mn, Nn from GreenTensor RCSCalculator (fresh instance, vacuum shell added)."""
    n_lay = len(radii) + 1
    calc = RCSCalculator(k0=k0, toch=n_terms + 1, n=n_lay,
                         a=list(radii) + [1.0],
                         eps=list(eps_arr) + [1.0],
                         miy=[1.0] * n_lay)
    calc.run_calculation()
    return calc.Mn[:n_terms], calc.Nn[:n_terms]


def sigma_r(x: float, a_n: np.ndarray, b_n: np.ndarray) -> float:
    """Normalized monostatic (backscatter) cross-section sigma_r / (pi a^2)."""
    n = np.arange(1, len(a_n) + 1)
    s = np.sum((2 * n + 1) * (-1.0) ** n * (a_n - b_n))
    return float(abs(s) ** 2) / x ** 2


def main() -> int:
    print("=" * 74)
    print("Example 6 — verification of the e^{+j omega t} time convention")
    print("=" * 74)
    failures = 0

    # --- V1: single layer -------------------------------------------------
    x, eps, N = 2 * math.pi, 2.5, 20
    a_c, b_c = mie_canon_multilayer(1.0, [x], [eps], [1.0], N)
    Mn, Nn = gt_coefficients(x, [1.0], [eps], N)
    e1 = max(float(np.max(np.abs(Mn - a_c))), float(np.max(np.abs(Nn - b_c))))
    ok = e1 < TOL
    failures += not ok
    print(f"V1 single layer (eps=2.5, k0a=2pi):  max|Mn-a_n^(h2)| = {e1:.3e}"
          f"   {'PASS' if ok else 'FAIL'}")

    # --- V2: 8-layer Fuchs Luneburg lens ----------------------------------
    k0, N2 = 4 * math.pi, 30
    a_c8, b_c8 = mie_canon_multilayer(k0, FUCHS_RADII, FUCHS_EPS, [1.0] * 8, N2)
    Mn8, Nn8 = gt_coefficients(k0, FUCHS_RADII, FUCHS_EPS, N2)
    e2 = max(float(np.max(np.abs(Mn8 - a_c8))), float(np.max(np.abs(Nn8 - b_c8))))
    ok = e2 < TOL
    failures += not ok
    print(f"V2 8-layer Fuchs LL (k0a=4pi):       max|Mn-a_n^(h2)| = {e2:.3e}"
          f"   {'PASS' if ok else 'FAIL'}")

    # --- V3: outgoing-wave direction --------------------------------------
    r = np.linspace(40.0, 41.0, 2001)      # far zone, k0 = 1
    zeta1 = r * (special.spherical_jn(1, r) - 1j * special.spherical_yn(1, r))
    dphase = np.gradient(np.unwrap(np.angle(zeta1)), r)
    slope = float(np.mean(dphase))          # expected ~ -k0 = -1
    ok = abs(slope + 1.0) < 1e-3
    failures += not ok
    print(f"V3 outgoing wave: d(arg zeta_1)/dr = {slope:+.6f} (expect -1.0)"
          f"      {'PASS' if ok else 'FAIL'}")
    print("   => with e^{+j omega t} the scattered phase fronts travel outward;"
          " h^{(1)} would give +1.0")

    # --- V4: energy conservation (lossless) --------------------------------
    uni = max(float(np.max(np.abs(Mn8.real - np.abs(Mn8) ** 2))),
              float(np.max(np.abs(Nn8.real - np.abs(Nn8) ** 2))))
    n = np.arange(1, N2 + 1)
    q_ext = (2 / k0 ** 2) * np.sum((2 * n + 1) * (Mn8.real + Nn8.real))
    q_sca = (2 / k0 ** 2) * np.sum((2 * n + 1) * (np.abs(Mn8) ** 2 + np.abs(Nn8) ** 2))
    e4 = abs(q_ext - q_sca) / q_ext
    ok = uni < 1e-10 and e4 < 1e-10
    failures += not ok
    print(f"V4 energy: max|Re Mn - |Mn|^2| = {uni:.3e};"
          f"  |Q_ext-Q_sca|/Q_ext = {e4:.3e}   {'PASS' if ok else 'FAIL'}")

    # --- V5: observable sweep ----------------------------------------------
    k0a_grid = np.linspace(1.0, 4 * math.pi, 60)
    sig_gt, sig_canon = [], []
    for x_i in k0a_grid:
        n_i = int(math.ceil(x_i + 4 * x_i ** (1 / 3) + 2)) + 4
        Mi, Ni = gt_coefficients(x_i, FUCHS_RADII, FUCHS_EPS, n_i)
        ac, bc = mie_canon_multilayer(x_i, FUCHS_RADII, FUCHS_EPS, [1.0] * 8, n_i)
        sig_gt.append(sigma_r(x_i, Mi, Ni))
        sig_canon.append(sigma_r(x_i, ac, bc))
    sig_gt, sig_canon = np.array(sig_gt), np.array(sig_canon)
    e5 = float(np.max(np.abs(sig_gt - sig_canon) / np.maximum(sig_canon, 1e-30)))
    ok = e5 < 1e-10
    failures += not ok
    print(f"V5 sigma_r(k0a) sweep, 60 pts:       max rel. diff = {e5:.3e}"
          f"   {'PASS' if ok else 'FAIL'}")

    # --- figure -------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    th = np.linspace(0, 2 * math.pi, 400)
    ax1.plot(0.5 + 0.5 * np.cos(th), 0.5 * np.sin(th), "k--", lw=0.8,
             label=r"unitarity $|w-\frac{1}{2}|=\frac{1}{2}$")
    ax1.plot(Mn8.real, Mn8.imag, "o", ms=6, mfc="none", color="tab:blue",
             label=r"GreenTensor $M_n$")
    ax1.plot(a_c8.real, a_c8.imag, "x", ms=5, color="tab:red",
             label=r"direct $h^{(2)}$: $a_n$")
    ax1.plot(Nn8.real, Nn8.imag, "s", ms=6, mfc="none", color="tab:green",
             label=r"GreenTensor $N_n$")
    ax1.plot(b_c8.real, b_c8.imag, "+", ms=7, color="tab:orange",
             label=r"direct $h^{(2)}$: $b_n$")
    ax1.set_xlabel(r"$\mathrm{Re}$"); ax1.set_ylabel(r"$\mathrm{Im}$")
    ax1.set_title(r"8-layer Fuchs LL, $k_0a=4\pi$: $M_n,N_n$ in $e^{+j\omega t}$ canon")
    ax1.axhline(0, color="gray", lw=0.5); ax1.axvline(0, color="gray", lw=0.5)
    ax1.set_aspect("equal"); ax1.legend(fontsize=8, loc="lower right")

    ax2.plot(k0a_grid, sig_gt, "-", color="tab:blue", lw=1.4,
             label="GreenTensor (hankel1 + conj)")
    ax2.plot(k0a_grid[::3], sig_canon[::3], "o", ms=4, mfc="none",
             color="tab:red", label=r"direct $h^{(2)}$ solver")
    ax2.set_xlabel(r"$k_0 a$")
    ax2.set_ylabel(r"$\sigma_r/(\pi a^2)$")
    ax2.set_title(r"Backscatter of the 8-layer Fuchs LL: identical observables")
    ax2.set_yscale("log"); ax2.grid(alpha=0.3); ax2.legend(fontsize=8)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(_HERE, f"time_convention_verification.{ext}"),
                    dpi=200, bbox_inches="tight")
    print(f"\nFigure saved: time_convention_verification.png / .pdf")

    print("=" * 74)
    print("ALL CHECKS PASS" if failures == 0 else f"{failures} CHECK(S) FAILED")
    print("=" * 74)
    return failures


if __name__ == "__main__":
    raise SystemExit(main())
