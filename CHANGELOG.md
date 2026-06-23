# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/) and the project adheres to
[Semantic Versioning](https://semver.org/).

## [0.2.0] — 2026-06-23

Major expansion of the original layered-sphere library into a full analytic
scattering suite, plus packaging cleanup for open-source release.

### Added
- **Unified public API** (`green_tensor.solvers`): one solver per geometry with a
  uniform interface — `SphereSolver`, `EllipsoidSolver`, `SpheroidSolver`,
  `CylinderSolver`, `ConeSolver`, and the `Cluster` (GMM) engine — plus thin
  `solve_*` facade functions. Bilingual (RU/EN) docstrings.
- **T-matrix interface** in the spherical VSWF basis (`tmatrix.py`,
  `DiagonalTMatrix`) shared by all geometries.
- **Vector spherical wave functions & translation** (`vswf.py`): Wigner-3j (Racah),
  Gaunt coefficients, VSWF M/N, and **closed-form analytic Cruzan
  translation-addition coefficients A, B** — removing the dense-packing limitation.
- **Generalized Multiparticle Mie engine** (`gmm.py`): self-consistent cluster
  solve, scattered field, far-field `C_sca`, optical-theorem `C_ext`.
- **Analytic geometry solvers**: triaxial ellipsoid and prolate/oblate spheroid
  (quasi-static depolarization, homogeneous + confocal-coated polarizability,
  Rayleigh cross sections, dipole T); infinite circular cylinder (exact 2-D Mie,
  TM/TE); cone via rigorous decomposition into a non-overlapping sphere cluster.
- **Geometry decomposition** (`decompose.py`): pack arbitrary bodies
  (sphere/box/cylinder/cone) into non-overlapping spheres for GMM.
- `Scatterer` protocol (`scatterer.py`) and `LayeredSphere` adapter over the core.
- Project documentation: MIT `LICENSE`, rewritten `README`, `CONTRIBUTING`,
  this changelog, and a real Sphinx site (replacing the placeholder).
- Independent-arbiter test suite covering every solver (`tests/run_all.py`).

### Changed
- Converted the package to proper relative intra-package imports; populated
  `green_tensor/__init__.py` with the public API and `__version__`.
- `pyproject.toml`: distribution `greentensor` ↔ import package `green_tensor`
  (`[tool.flit.module]`), MIT license metadata, `docs`/`dev` extras.

### Moved
- Original standalone research scripts (`02_cilindre.py`, `03_decart_*.py`,
  `Bistatic_RCS_*.py`, `lin_polar.py`, `examples.json`) → `green_tensor/legacy/`
  (archival, not part of the public API).

### Removed
- Stray `greentensor.py` version stub (version now lives in the package).
- Generated `outputs/` artifacts from version control (now git-ignored).

### Preserved
- `green_tensor/01_sphere.py` — the canonical exact sphere core (TGF/Mie). Its
  mathematics is unchanged.

## [0.1.0]

- Initial open-source release: layered-sphere electromagnetic scattering by the
  tensor Green's function method; Luneburg/Maxwell-lens examples; basis of the two
  IEEE publications.
