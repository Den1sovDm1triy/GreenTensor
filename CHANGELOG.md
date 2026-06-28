<!-- SPDX-License-Identifier: MIT -->

# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/) and the project adheres to
[Semantic Versioning](https://semver.org/).

## [0.5.0] — 2026-06-28

Сокращение области охвата до **только точных аналитических решений метода тензорных
функций Грина (ТФГ)**. Цель — научная библиотека, в которой каждый решатель строго
аналитичен и верифицируем; приближённые и чисто численные методы исключены.

### Removed (BREAKING)
- **EBCM / метод нулевого поля** и несферические примитивы `Spheroid`,
  `FiniteCylinder`, `Cone` (модуль `green_tensor/ebcm.py`) — численный метод,
  чувствительный к рёбрам/вершинам и неустойчивый для слоистых тел.
- **Разложение тела в кластер сфер** (`green_tensor/decompose.py`) и решатели
  `ConeSolver`, `FiniteCylinderSolver` — приближение (упаковка сферами +
  поправка Максвелла–Гарнетта).
- **Квазистатика Рэлея**: `green_tensor/ellipsoid.py`, `green_tensor/spheroid.py`,
  класс `EllipsoidSolver` и функция `solve_ellipsoid` — точны лишь в пределе
  `k·a ≪ 1`, не полноволновые.
- **`green_tensor/mie_core.py`** — дубликат сферы; единственный путь сферы теперь —
  каноническое ядро `01_sphere.py` через фасад `sphere_core`.
- **Вращения Вигнера T-матрицы** (`vswf.rotate_tmatrix`, `wigner_d`, `wigner_D`,
  `rotation_matrix`) — служили только ориентации EBCM-тел; сферы изотропны, кластеру
  нужна лишь трансляция.
- Удалены связанные тесты и примеры (`composite_primitives_demo.py`,
  `complex_geometry_demo.py`).

### Kept (точное аналитическое ТФГ)
- **Слоистая сфера** — точное ядро `01_sphere.py` (`SphereSolver`, `solve_sphere`).
- **Бесконечный цилиндр** — точная 2D-аналитика, одно-/многослойный, нормальное и
  косое падение (`CylinderSolver`, `LayeredCylinderSolver`, ТФГ/эквивалентные линии).
- **Кластер невзаимопересекающихся сфер** — строгая суперпозиция через теорему
  сложения Крузана–Стейна (`Cluster`, GMM).

### Added
- **Антенный режим сферы** (источник/рупор на поверхности) реализован прямо в каноне
  `01_sphere.py` — переключатель `problem="diffraction"|"antenna"` (ранее формула рупора
  была в коде закомментирована). Проброшен через `sphere_core`,
  `SphereSolver.pattern(..., problem="antenna")`, webapp и Studio. Это задача излучения:
  выводится диаграмма направленности (сечения рассеяния и ЭПР для неё не определены).

### Changed
- Публичный API сведён к `SphereSolver`, `CylinderSolver`, `LayeredCylinderSolver`,
  `Cluster` + обёртки `solve_sphere`/`solve_cylinder`/`solve_layered_cylinder`/`solve_cluster`.
- Документация и `GreenTensor_Theory.tex` приведены к новому охвату (удалены разделы
  EBCM, конечного цилиндра, конуса, сфероида SVM, триаксиального эллипсоида,
  разложения в сферы). Версия `0.5.0`.

## [0.4.0] — 2026-06-27

### Added
- **Rigorous full-wave EBCM / null-field T-matrix solver** for non-spherical
  axisymmetric bodies (`green_tensor/ebcm.py`): spheroid, finite cylinder, cone,
  including **layered** bodies via the Peterson–Ström recursion. Calibrated to the
  layered-sphere core (`sphere_core`) and validated for energy conservation, the
  Rayleigh limit, and reduction to Mie.
- **EBCM primitives exposed in the public API**: `green_tensor.Spheroid`,
  `green_tensor.FiniteCylinder`, `green_tensor.Cone` — implement the `Scatterer`
  protocol and assemble through `Cluster` (GMM) for exact single-body or cluster
  cross sections. This is now the **main route for non-spherical scattering**.
- **Arbitrary orientation** of EBCM primitives via Wigner-D rotation of the
  T-matrix (`vswf.rotate_tmatrix`, `euler=(α,β,γ)` on each primitive), validated by
  cross-section covariance.
- GMM generalized to **full (non-diagonal) T-matrices** (`gmm._full_t`); spheres
  remain bit-for-bit unchanged.

### Changed (BREAKING)
- **Non-spherical solvers consolidated onto the rigorous EBCM/TGF family.** The
  approximate quasi-static `SpheroidSolver` and `solve_spheroid` are **removed**:
  the full-wave spheroid is now `green_tensor.Spheroid` (EBCM), and the Rayleigh
  limit is available via `EllipsoidSolver(a_eq, a_eq, c_ax, …)` (spheroid = ellipsoid
  with b = a — numerically identical depolarization/polarizability).
- `green_tensor/cone.py` **removed**: its geometry helpers (`cone_indicator`,
  `decompose_cone`) moved into `green_tensor/decompose.py`; `ConeSolver` now wires
  the decompose-fallback path, while the rigorous cone is `green_tensor.Cone` (EBCM).
- Retained for their **unique physics** (not covered by EBCM): `EllipsoidSolver`
  (triaxial, quasi-static), `CylinderSolver`/`LayeredCylinderSolver` (infinite,
  layered, oblique 2D), `decompose` (arbitrary-blob fallback), `mie_core`
  (antenna-mode + verified oracle).

### Fixed
- `mie_core` magneto-dielectric (μ ≠ 1) interface factors — magnetic sphere now
  matches the Kerker closed-form arbiter; GMM axial-incidence NaN at the poles.

### Removed
- Out-of-scope `examples/Example 1 …/test.py` (Yandex-GPT HTML analyzer, not
  electrodynamics) — its previously committed API key must be revoked in the
  Yandex Cloud console.

## [0.3.0] — 2026-06-23

### Added
- **Layered-cylinder scattering solver** via the tensor Green's function (TGF) /
  equivalent transmission-line method (Daylis & Shabunin 2017):
  - `cylinder.layered_coeff` / `cross_sections_layered` — radially multilayer
    cylinder at **normal incidence** (decoupled TM/TE transfer-matrix recursion).
  - `cylinder.layered_coeff_oblique` / `cross_sections_oblique` — **oblique
    incidence** (h ≠ 0) with E–H polarization coupling via a 4×4 layer transfer
    matrix; returns co- and cross-polarized scattering.
  - `LayeredCylinderSolver` (with `theta` for oblique) + `solve_layered_cylinder`
    in the public API.
  - `tests/analytic_cylinder.py`: independent arbiters — Bohren–Huffman (normal)
    and Wait (1955) / Kavaklıoğlu (oblique). New solver verified to machine
    precision against both, plus subdivision invariance and energy conservation.
  - Theory: new section «Слоистый бесконечный цилиндр (ТФГ)» in
    `GreenTensor_Theory.tex` (+ Daylis–Shabunin, Wait, Kavaklıoğlu references).
- **Finite cylinder integrated into the GMM engine via decomposition**:
  `decompose.decompose_cylinder` packs a finite cylinder into a non-overlapping
  sphere cluster; `FiniteCylinderSolver` (decompose → `Cluster`) added to the public
  API, mirroring `ConeSolver`. Verified non-overlapping/inside + GMM integration.
- **Metal guard for the decomposition route** (`decompose.is_metal_layer` /
  `reject_metal_packing`): sphere packing of a metal body creates spurious
  inter-sphere cavity re-reflections (a metallic "sponge"), so `decompose_cylinder`,
  `decompose_cone`, and the `FiniteCylinderSolver`/`ConeSolver` `.decompose` methods
  now raise `ValueError` when the **outermost** layer is metallic (`Re(eps) < 0` or
  skin depth `< r_sphere`, via an optional `k`). A metal core enclosed by a dielectric
  outer layer is allowed; `allow_metal=True` overrides for a genuine metal-sphere cluster.
- **Decomposition optimizer — Maxwell–Garnett effective-medium correction**: packing
  leaves vacuum gaps that lower the cluster's effective permittivity; `effective_medium=True`
  on the `decompose` route chooses the sphere ε (closed-form `maxwell_garnett_eps`) so the
  effective medium of (spheres at filling `f` in vacuum) equals the body's true ε. Exact in
  the quasi-static limit (`Σα_corrected == α_solid`; uncorrected is biased by `f`).
  Filling-limited (max reachable ε = `(1+2f)/(1−f)`; raises if `f` too low — denser packing
  needed). Added `maxwell_garnett_effective`, `maxwell_garnett_eps`, `effective_medium_eps`,
  and `packing_report` (filling, overlap margin, GMM size estimate).
- **FCC packing lattice** (`pack_spheres(..., lattice="fcc")`): face-centered-cubic gives a
  ~√2× denser, still non-overlapping packing (higher filling, better boundary conformity),
  raising the Maxwell–Garnett reachable ε. Threaded through `decompose_cylinder`,
  `decompose_cone`, and the `FiniteCylinderSolver`/`ConeSolver` `.decompose` methods
  (`lattice="cubic"` default keeps prior behavior).

### Fixed
- `cylinder._coeffs` had TM/TE labels swapped vs the Bohren–Huffman convention
  (TM = E parallel to the axis); corrected. Cross-section/energy results unaffected.

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
