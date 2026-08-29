# SPDX-License-Identifier: MIT
"""Example 10 — neural-network identification of a coating from its spectrum.

Demonstrates how GreenTensor turns the *inverse* radar-diagnostics task into a
supervised-learning problem: the library computes back-scattering spectra of a
PEC sphere with a magnetodielectric coating (the forward map), and a neural
network learns to invert them — recovering coating thickness and complex
(eps, mu) from the spectrum.

This script builds a tiny labelled sample (a few configurations) to show the
forward map that the network inverts. Strongly absorbing coatings make
GreenTensor switch to its arbitrary-precision branch, so a handful of spectra
is enough for a demo. The full, reproducible dataset (6000 configurations,
generated with a fast bundled solver), the reference inverse network and the
robustness study are published separately:

    https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset

Run:  python3 nn_coating_identification.py
"""
from __future__ import annotations
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from green_tensor.multilayer_mie import mie_pec_with_coating
from green_tensor.mie_reference import cross_sections, wiscombe_nmax

# 32-point spectrum over the resonance region (as in the published dataset)
GRID = np.linspace(1.5, 8.5, 32)          # k0 * R

# A few representative coatings (moderate losses -> fast double-precision path).
# Each row is (t/R, eps', eps'', mu', mu''). The published dataset spans the full
# parameter box, including strongly absorbing coatings, via a fast bundled solver.
CASES = [
    ("thin dielectric",        [0.04,  3.4, 0.03, 1.0, 0.0]),
    ("medium dielectric",      [0.08,  6.0, 0.20, 1.0, 0.0]),
    ("mild magnetodielectric", [0.05,  8.0, 0.50, 1.5, 0.4]),
    ("thicker mild-loss",      [0.10,  9.0, 0.50, 2.0, 0.4]),
]


def backscatter_spectrum(p):
    """Forward map (GreenTensor): coating parameters -> lg Q_back spectrum.

    p = (t/R, eps', eps'', mu', mu'');  engineering convention e^{+jwt}:
    eps = eps' - j eps'',  mu = mu' - j mu''.
    """
    t, eps, mu = p[0], p[1] - 1j * p[2], p[3] - 1j * p[4]
    r_out = 1.0 + t
    n_max = wiscombe_nmax(abs(np.sqrt(eps * mu)) * GRID[-1] * r_out) + 12
    s = np.empty(len(GRID))
    for i, k0r in enumerate(GRID):
        a, b = mie_pec_with_coating(k0r, 1.0, [r_out], [eps], [mu], n_max)
        s[i] = cross_sections(a, b, k0r * r_out)[2]      # Q_back
    return np.log10(np.maximum(s, 1e-12))


def main():
    print("Forward map computed by GreenTensor (coating parameters -> spectrum):\n")
    print(f"{'coating':<24}{'lg Q_back min':>14}{'lg Q_back max':>14}")
    S, P = [], []
    for name, p in CASES:
        s = backscatter_spectrum(p)
        S.append(s); P.append(p)
        print(f"{name:<24}{s.min():>14.2f}{s.max():>14.2f}")
    S, P = np.array(S), np.array(P)

    print("\nThe inverse problem the neural network solves:")
    print("  given the 32-point spectrum S  ->  recover the 5 parameters P")
    print("  (t/R, eps', eps'', mu', mu').")
    print("\nFull dataset (6000 configurations), the reference inverse MLP")
    print("(R2 ~ 0.98 for thickness and magnetic parameters; the dielectric")
    print("loss eps'' is the hard one) and the robustness study are here:")
    print("  https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset")
    return S, P


if __name__ == "__main__":
    main()
