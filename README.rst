.. SPDX-License-Identifier: MIT

GreenTensor
===========

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: LICENSE
   :alt: MIT License

**GreenTensor** is a Python library for analytic electromagnetic-scattering
analysis on heterogeneous structures. Its scope is restricted to **exact analytic**
tensor Green's function (TGF) solutions: the *exact* radially layered sphere
(TGF / Mie), the infinite layered/oblique circular cylinder (exact 2-D analytics),
and their rigorous superposition into a **cluster of non-overlapping spheres** via
the **Generalized Multiparticle Mie (GMM)** method (Cruzan–Stein addition theorem),
sharing one **T-matrix interface** in the spherical vector-wave-function (VSWF)
basis. Approximate or purely numerical methods are intentionally excluded.

The library underpins two IEEE publications:

- `GreenTensor: Efficient Electromagnetic Scattering Analysis on Heterogeneous Spherical Structures <https://ieeexplore.ieee.org/document/10615009>`_
- `GreenTensor: A Python Library for Electromagnetic Wave Propagation Analysis <https://ieeexplore.ieee.org/document/10584060>`_

The complete mathematical derivation (per geometry, verified against the
literature) is in ``GreenTensor_Theory.tex`` at the repository root.

Features
--------

- **Exact layered-sphere core** (Mie / TGF): multilayer, complex permittivity
  (absorption), magnetic layers, PEC limit, linear and circular polarization, both
  the plane-wave **diffraction** and surface-source (**antenna**) problems.
- **Exact analytic solvers only:**

  - ``SphereSolver`` — layered sphere, exact (Mie / TGF).
  - ``CylinderSolver`` / ``LayeredCylinderSolver`` — infinite cylinder, exact 2-D
    (homogeneous, layered, normal + oblique incidence; TGF / transmission-line).
  - ``Cluster`` — self-consistent GMM assembly of **non-overlapping spheres**
    (Cruzan–Stein addition theorem; closed-form coefficients valid at any
    separation).

- **Unified T-matrix interface** (spherical VSWF) feeding the GMM engine.
- **Honest scope.** Only exact, verifiable analytics are shipped. Non-separable
  geometries with no exact closed-form solution (e.g. the finite cylinder) raise
  ``NotImplementedError`` rather than returning unverified results; approximate or
  numerical methods (EBCM for edged bodies, Rayleigh quasi-statics, sphere-packing
  decomposition) are not included.

Open datasets
-------------

Two open benchmark datasets are generated with this library — the same physics
taken in both directions:

- **LL-Surrogate-6000** — forward (surrogate) map of a four-layer Luneburg lens:
  ``k0a, eps1..eps4`` → ``Q_sca``, ``log10 Q_back``; 6000 rows, labels from
  ``green_tensor.sphere_core.MieSphere``.
  `luneburg-lens-nn-dataset <https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset>`_
- **RPM-Inverse-6000** — inverse identification of a magnetodielectric coating on
  a PEC sphere: a 32-point back-scattering spectrum → ``t/R``, ``eps'``, ``eps''``,
  ``mu'``, ``mu''``; 6000 rows, self-contained solver bundled.
  `rpm-coating-inverse-dataset <https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset>`_

.. figure:: https://raw.githubusercontent.com/Den1sovDm1triy/luneburg-lens-nn-dataset/main/figures/gt-surrogate-parity.png
   :width: 100%
   :alt: Surrogate predictions against GreenTensor labels
   :target: https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset

   LL-Surrogate-6000 — a plain MLP against the GreenTensor labels: ``R² = 0.992``
   for ``Q_sca``, ``0.804`` for the deliberately sharp ``log10 Q_back``.

.. figure:: https://raw.githubusercontent.com/Den1sovDm1triy/rpm-coating-inverse-dataset/main/figures/spectra-examples.png
   :width: 100%
   :alt: Back-scattering spectra of a bare PEC sphere and three coatings
   :target: https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset

   RPM-Inverse-6000 — the input signal: back-scattering of a bare PEC sphere
   against three coating types; the network reads exactly this reshaping.

Both are MIT-licensed, fixed-seed reproducible and ship reference models,
robustness scripts and figures. Overview, schemas and verification recipe:
`DATASETS.md <DATASETS.md>`_; a runnable tour is
``examples/Example 12 - Open Datasets``.

Installation
------------

.. code-block:: bash

   pip install .            # from a checkout
   # optional extras:
   pip install ".[docs]"    # Sphinx docs
   pip install ".[dev]"     # pytest + coverage

Requires Python ≥ 3.9, NumPy, SciPy, Matplotlib.

The distribution is named ``greentensor``; the import package is ``green_tensor``.

Quick start
-----------

.. code-block:: python

   import green_tensor as gt

   # Layered sphere — exact cross sections (normalized to πR²)
   gt.solve_sphere(radius=1.0, eps=[2.25], k=3.0)
   # -> {'q_sca': ..., 'q_ext': ..., 'q_abs': ..., 'q_back': ...}

   # Object API: T-matrix, scattering pattern
   sphere = gt.SphereSolver(radius=1.0, eps=[2.25])
   T = sphere.t_matrix(k=3.0)            # DiagonalTMatrix in the VSWF basis
   diagram = sphere.pattern(k=3.0)       # scattering diagram

   # Infinite layered cylinder (exact 2-D, TGF) — normal & oblique incidence
   gt.solve_layered_cylinder(radius=1.0, eps=[4 + 0.5j, 2.0], a_norm=[0.6, 1.0], k=2.0)

   # Cluster of non-overlapping spheres (GMM)
   s1 = gt.SphereSolver(0.5, [2.25], position=(0, 0, -1.5)).as_scatterer()
   s2 = gt.SphereSolver(0.5, [2.25], position=(0, 0, +1.5)).as_scatterer()
   gt.Cluster([s1, s2]).cross_sections(k=3.0, khat=(1, 0, 0), pol=(0, 0, 1), nmax=6)

See ``examples/`` for end-to-end scripts (Luneburg / Maxwell lenses).

Convention
----------

Time dependence ``e^{-iωt}``; outgoing waves ``~ h_n^{(1)}``; vacuum host medium;
``k`` is the host wavenumber and the size parameter is ``x = k · radius``.
``Im(eps) > 0`` denotes absorption.

Project layout
--------------

.. code-block:: text

   green_tensor/        importable package
     01_sphere.py       canonical exact sphere core (TGF/Mie) — math preserved
     sphere_core.py     importable facade over 01_sphere.py
     tmatrix.py         diagonal T-matrix + cross sections
     vswf.py            VSWF, Wigner-3j/Gaunt, closed-form translation-addition
     scatterer.py       Scatterer protocol + LayeredSphere
     gmm.py             Generalized Multiparticle Mie assembly engine (sphere clusters)
     cylinder.py        infinite/layered/oblique cylinder (exact 2-D, TGF)
     solvers.py         unified public API (SphereSolver, CylinderSolver,
                        LayeredCylinderSolver, Cluster, solve_*)
     legacy/            archival original research scripts (not maintained)
   tests/               independent-arbiter test suite (run_all.py)
   docs/                Sphinx documentation
   examples/            worked examples
   webapp/              browser front-end for the layered-sphere core
   GreenTensor_Theory.tex   full mathematical derivation

Testing
-------

.. code-block:: bash

   python3 tests/run_all.py      # no pytest required
   # or
   pytest tests/

Tests validate every solver against independent analytic references
(Bohren & Huffman Mie arbiter, energy conservation, optical theorem).

Documentation
-------------

Built with Sphinx (``docs/``); hosted on `Read the Docs
<https://greentensor.readthedocs.io>`_.

.. code-block:: bash

   pip install ".[docs]"
   sphinx-build -b html docs docs/_build/html

License
-------

MIT — see `LICENSE <LICENSE>`_. © 2025–2026 D.V. Denisov, V.Ya. Noskov,
I.O. Skumatenko (Ural Federal University). Every source file carries an
``SPDX-License-Identifier: MIT`` header.

Citation
--------

If you use GreenTensor in academic work, please cite both the software and the
underlying IEEE publications. Machine-readable metadata lives in
`CITATION.cff <CITATION.cff>`_ (GitHub's *“Cite this repository”* button).

Software (BibTeX):

.. code-block:: bibtex

   @software{greentensor,
     title   = {{GreenTensor}: analytic electromagnetic scattering on a layered-sphere core},
     author  = {Denisov, Dmitriy V. and Noskov, Vitaliy Ya. and Skumatenko, Ilya O.},
     version = {0.5.0},
     year    = {2026},
     license = {MIT},
     url     = {https://github.com/Den1sovDm1triy/GreenTensor},
   }

Underlying methodology — the two IEEE papers linked at the top of this README.
