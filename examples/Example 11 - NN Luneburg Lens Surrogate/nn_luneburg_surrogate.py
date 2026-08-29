# SPDX-License-Identifier: MIT
"""Example 11 — neural-network surrogate of a multilayer Luneburg lens.

Demonstrates the *forward* (surrogate) direction: GreenTensor computes the
scattering cross sections of a four-layer Luneburg lens, and a neural network
learns to reproduce them directly from the layer permittivities — replacing the
full-wave solver inside an optimization loop by a microsecond lookup.

    parameters (k0a, eps1..eps4)  ->  (Q_sca, log10 Q_back)

This script evaluates a few lens configurations to show the map the network
approximates. The full, reproducible dataset (6000 configurations), the
reference surrogate network and the scripts are published separately:

    https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset

Run:  python3 nn_luneburg_surrogate.py
"""
from __future__ import annotations
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from green_tensor.sphere_core import MieSphere

RADII = [0.25, 0.5, 0.75, 1.0]            # equal-step normalized layer boundaries

# A few four-layer Luneburg-type lenses: eps decreasing from centre to surface.
# Row = (k0a, eps1, eps2, eps3, eps4).
CASES = [
    ("ideal Luneburg, k0a=4pi", [4 * np.pi, 1.94, 1.66, 1.31, 1.06]),
    ("ideal Luneburg, k0a=6pi", [6 * np.pi, 1.94, 1.66, 1.31, 1.06]),
    ("flatter profile",         [4 * np.pi, 1.70, 1.55, 1.35, 1.15]),
    ("steeper profile",         [5 * np.pi, 2.00, 1.60, 1.25, 1.05]),
]


def cross_sections(p):
    """Forward map (GreenTensor): (k0a, eps1..eps4) -> (Q_sca, log10 Q_back)."""
    k0a, eps = p[0], list(p[1:])
    cs = MieSphere(k0=k0a, a=RADII, eps=eps, problem="diffraction").cross_sections()
    return cs["q_sca"], np.log10(max(cs["q_back"], 1e-8))


def main():
    print("Forward map computed by GreenTensor (lens profile -> cross sections):\n")
    print(f"{'configuration':<26}{'Q_sca':>10}{'log10 Q_back':>14}")
    X, Y = [], []
    for name, p in CASES:
        qs, lqb = cross_sections(p)
        X.append(p); Y.append([qs, lqb])
        print(f"{name:<26}{qs:>10.3f}{lqb:>14.3f}")
    X, Y = np.array(X), np.array(Y)

    print("\nThe surrogate the neural network learns:")
    print("  given (k0a, eps1..eps4)  ->  predict (Q_sca, log10 Q_back).")
    print("\nFull dataset (6000 configurations), the reference surrogate MLP")
    print("(R2 ~ 0.99 for Q_sca) and the scripts are here:")
    print("  https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset")
    return X, Y


if __name__ == "__main__":
    main()
