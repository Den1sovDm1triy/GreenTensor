# SPDX-License-Identifier: MIT
"""Example 12 — a tour of the two open datasets built with GreenTensor.

The library ships two benchmark datasets as separate repositories:

    LL-Surrogate-6000   forward  (k0a, eps1..eps4) -> (Q_sca, log10 Q_back)
                        https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset
    RPM-Inverse-6000    inverse  32-point lg Q_back spectrum -> (t/R, eps', eps'', mu', mu'')
                        https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset

This script finds them on disk, prints what is inside, and — the point of the
example — *verifies* the stored labels by recomputing them: the Luneburg set
against this library, the coating set against the solver bundled with it.

Where the datasets are looked for, in order:

    1. $GT_LL_DATASET / $GT_RPM_DATASET          (explicit paths)
    2. ./luneburg-lens-nn-dataset, ./rpm-coating-inverse-dataset   (cwd)
    3. the same names next to this repository and next to this example

Nothing is downloaded: if a dataset is missing, the script says how to clone it
and moves on. See DATASETS.md at the repository root for the overview.

Run:  python3 open_datasets_tour.py
"""
from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, REPO)

LL_NAME = "luneburg-lens-nn-dataset"
RPM_NAME = "rpm-coating-inverse-dataset"
LL_URL = f"https://github.com/Den1sovDm1triy/{LL_NAME}"
RPM_URL = f"https://github.com/Den1sovDm1triy/{RPM_NAME}"

# Fixed geometry of the Luneburg set: normalized outer radii of the four layers.
LL_RADII = [0.25, 0.5, 0.75, 1.0]


def find_dataset(name: str, env_var: str, marker: str) -> str | None:
    """First directory that holds `marker`, searched in the documented order."""
    candidates = [os.environ.get(env_var), os.path.join(os.getcwd(), name),
                  os.path.join(REPO, name), os.path.join(os.path.dirname(REPO), name),
                  os.path.join(HERE, name)]
    for path in candidates:
        if path and os.path.isfile(os.path.join(path, marker)):
            return os.path.abspath(path)
    return None


def missing(name: str, url: str, env_var: str) -> None:
    print(f"  not found. Get it with:\n    git clone {url}")
    print(f"  or point the script at an existing copy:  export {env_var}=/path/to/{name}")


def tour_luneburg() -> None:
    """Forward set: describe it, then recompute a few labels with the library."""
    print("\n" + "=" * 72)
    print("LL-Surrogate-6000 — forward map of a four-layer Luneburg lens")
    print("=" * 72)

    path = find_dataset(LL_NAME, "GT_LL_DATASET", "ll_surrogate_6000.npz")
    if path is None:
        missing(LL_NAME, LL_URL, "GT_LL_DATASET")
        return
    print(f"  dataset: {path}")

    d = np.load(os.path.join(path, "ll_surrogate_6000.npz"))
    X, Y = d["X"], d["Y"]
    print(f"  X {X.shape}  (k0a, eps1..eps4)      k0a in "
          f"[{X[:, 0].min():.2f}, {X[:, 0].max():.2f}]")
    print(f"  Y {Y.shape}  (Q_sca, log10 Q_back)  Q_sca in "
          f"[{Y[:, 0].min():.3f}, {Y[:, 0].max():.3f}], "
          f"log10 Q_back in [{Y[:, 1].min():.2f}, {Y[:, 1].max():.2f}]")

    try:
        from green_tensor.sphere_core import MieSphere
    except Exception as exc:                      # library not importable
        print(f"  (skipping the check: GreenTensor is not importable — {exc})")
        return

    print("\n  Verification — recompute the stored labels with GreenTensor:")
    print(f"    {'row':>6}{'Q_sca (data)':>15}{'Q_sca (solver)':>16}"
          f"{'log10 Q_back':>15}{'(solver)':>12}")
    worst = 0.0
    for i in (0, 1, 2, 1234, 5999):
        cs = MieSphere(k0=X[i, 0], a=LL_RADII, eps=list(X[i, 1:]),
                       problem="diffraction").cross_sections()
        q_sca, lq_back = cs["q_sca"], np.log10(max(cs["q_back"], 1e-8))
        worst = max(worst, abs(q_sca - Y[i, 0]), abs(lq_back - Y[i, 1]))
        print(f"    {i:>6}{Y[i, 0]:>15.4f}{q_sca:>16.4f}{Y[i, 1]:>15.4f}{lq_back:>12.4f}")
    print(f"  largest deviation on these rows: {worst:.2e}  "
          f"({'labels reproduce exactly' if worst < 1e-9 else 'check the solver version'})")

    print("\n  What to do with it:")
    print("    - train your own surrogate on X -> Y; the reference NumPy MLP is")
    print(f"      `python3 train_baseline.py 6000` inside {LL_NAME}")
    print("    - benchmark against R2 ~ 0.99 (Q_sca) and ~ 0.80 (log10 Q_back)")
    print("    - the sharp target is deliberately hard: frequency encoding or")
    print("      physics-informed regularization are the natural next steps")


def tour_coating() -> None:
    """Inverse set: describe it, then recompute a spectrum with its own solver."""
    print("\n" + "=" * 72)
    print("RPM-Inverse-6000 — coating identification from a back-scattering spectrum")
    print("=" * 72)

    path = find_dataset(RPM_NAME, "GT_RPM_DATASET", "rpm_inverse_6000.npz")
    if path is None:
        missing(RPM_NAME, RPM_URL, "GT_RPM_DATASET")
        return
    print(f"  dataset: {path}")

    d = np.load(os.path.join(path, "rpm_inverse_6000.npz"))
    S, P, grid = d["S"], d["P"], d["GRID"]
    print(f"  S {S.shape}  lg Q_back on a {len(grid)}-point grid, "
          f"k0R in [{grid[0]:.1f}, {grid[-1]:.1f}]")
    print(f"  P {P.shape}  (t/R, eps', eps'', mu', mu'')")
    print(f"    t/R    in [{P[:, 0].min():.3f}, {P[:, 0].max():.3f}]")
    print(f"    eps'   in [{P[:, 1].min():.2f}, {P[:, 1].max():.2f}], "
          f"eps'' in [{P[:, 2].min():.2f}, {P[:, 2].max():.2f}]")
    print(f"    mu'    in [{P[:, 3].min():.2f}, {P[:, 3].max():.2f}], "
          f"mu''  in [{P[:, 4].min():.2f}, {P[:, 4].max():.2f}]")

    solver_file = os.path.join(path, "forward_solver.py")
    if not os.path.isfile(solver_file):
        print("  (skipping the check: forward_solver.py is not in the dataset copy)")
        return

    spec = importlib.util.spec_from_file_location("rpm_forward_solver", solver_file)
    solver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(solver)

    print("\n  Verification — recompute stored spectra with the bundled solver:")
    worst = 0.0
    for i in (0, 5, 4321):
        t, eps_re, eps_im, mu_re, mu_im = P[i]
        spectrum, _ = solver.backscatter_spectrum(
            grid, t, complex(eps_re, -eps_im), complex(mu_re, -mu_im))
        dev = float(np.max(np.abs(np.asarray(spectrum) - S[i])))
        worst = max(worst, dev)
        print(f"    row {i:>5}: max |recomputed - stored| = {dev:.2e}")
    print(f"  largest deviation: {worst:.2e}  "
          f"({'spectra reproduce to machine precision' if worst < 1e-9 else 'check the solver version'})")

    print("\n  What to do with it:")
    print("    - train the inverse map S -> P; the reference MLP is")
    print(f"      `python3 train_inverse.py` inside {RPM_NAME}")
    print("    - it also reports the metrics under 0.5 dB spectrum noise —")
    print("      that noise test is the honest measure for a diagnostics model")
    print("    - `robustness_study.py` gives 5-fold CV, seed spread, learning curve")


def main() -> int:
    print("Open datasets built with GreenTensor — see DATASETS.md at the repo root.")
    tour_luneburg()
    tour_coating()
    print("\nBoth datasets are MIT-licensed, fixed-seed reproducible, and can be")
    print("regenerated from scratch: `generate_dataset.py` in each repository.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
