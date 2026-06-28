<!-- SPDX-License-Identifier: MIT -->

# Contributing to GreenTensor

Thanks for your interest in GreenTensor! This guide covers the development setup,
conventions, and the one hard rule about the mathematical core.

## Development setup

```bash
git clone https://github.com/Den1sovDm1triy/GreenTensor
cd GreenTensor
pip install -e ".[dev,docs]"
```

The distribution is `greentensor`; the import package is `green_tensor`.

## Running the tests

The suite needs **no pytest** — it ships a self-contained runner with an
independent analytic arbiter (Bohren & Huffman Mie):

```bash
python3 tests/run_all.py     # all checks
pytest tests/                # equivalent, if pytest is installed
```

Every change must keep the full suite green. New physics must be validated
against an **independent** reference (closed-form limit, energy conservation,
optical theorem, or a published benchmark) — not against the implementation
itself.

## The one hard rule: do not change the sphere core math

`green_tensor/01_sphere.py` is the **canonical mathematical core** (tensor
Green's function / Mie recursion). Its mathematical lines are frozen: the
verified correctness fixes (external-medium closure, impedance recursion range,
complex refractive index, metal-phase handling, scaled Bessel functions) are part
of the canon. **Do not refactor, "simplify", or migrate that math.** The
importable mirror is `green_tensor/mie_core.py`; keep the two consistent.

If you believe the core math needs to change, open an issue first and get
explicit sign-off before touching those lines.

## No fake or unverified code

Branches that would require special-function machinery not available in SciPy
(full-wave confocal spheroid, finite cylinder, sphero-conal cone) raise a
documented `NotImplementedError` pointing at the relevant theory section. Do not
replace them with placeholder/approximate code presented as exact. Complex
geometry is handled rigorously via `decompose` → `Cluster` (GMM) over the
analytic sphere family.

## Adding a new solver

1. Implement the geometry's analytics in its own module (mirror the style of
   `ellipsoid.py` / `cylinder.py`), keeping functions small and pure.
2. To make it composable in clusters, emit a T-matrix in the spherical VSWF basis
   and/or expose a `Scatterer`-protocol object (see `scatterer.py`).
3. Add a uniform wrapper to `green_tensor/solvers.py` (an OO `*Solver` class plus,
   optionally, a `solve_*` facade function) with **bilingual** docstrings.
4. Re-export it from `green_tensor/__init__.py` (`__all__`).
5. Add `tests/test_<geometry>.py` validating against an independent reference,
   plus a `tests/test_solvers.py` case asserting the wrapper reproduces the
   underlying verified function exactly. Wire both into `tests/run_all.py`.
6. Document it in `docs/solvers.rst` and `docs/api.rst`.

## Coding style

- Many small, focused files (≤ ~400 lines); high cohesion, low coupling.
- Prefer immutable, pure functions; validate inputs at boundaries; never swallow
  errors silently.
- Intra-package imports are **relative** (`from . import vswf`).
- Public API docstrings are **bilingual** (RU + EN); inline math docstrings may
  stay in Russian to match the existing core.
- Follow the conventions in `GreenTensor_Theory.tex` (`e^{-iωt}`, outgoing
  `h_n^{(1)}`, vacuum host, size parameter `x = k·radius`).

## Commits & PRs

- Conventional commits (`feat:`, `fix:`, `docs:`, `test:`, `refactor:`, `chore:`).
- Keep the suite green and update docs/CHANGELOG with user-facing changes.
