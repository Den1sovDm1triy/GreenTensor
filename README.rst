.. SPDX-License-Identifier: MIT
.. Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

GreenTensor
===========

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: LICENSE
   :alt: MIT License

**GreenTensor** is a Python library for analytic electromagnetic-scattering
analysis on heterogeneous structures. Its mathematical core is the *exact*
solution for a radially layered sphere obtained by the **tensor Green's function
(TGF / Mie) method**. On top of that core, a family of analytic solvers
(spheroid, ellipsoid, cylinder, cone) shares a single **T-matrix interface** in
the spherical vector-wave-function (VSWF) basis, and a single **Generalized
Multiparticle Mie (GMM)** engine composes them into solutions for arbitrarily
complex geometry.

The library underpins two IEEE publications:

- `GreenTensor: Efficient Electromagnetic Scattering Analysis on Heterogeneous Spherical Structures <https://ieeexplore.ieee.org/document/10615009>`_
- `GreenTensor: A Python Library for Electromagnetic Wave Propagation Analysis <https://ieeexplore.ieee.org/document/10584060>`_

The complete mathematical derivation (per geometry, verified against the
literature) is in ``GreenTensor_Theory.tex`` at the repository root.

Features
--------

- **Exact layered-sphere core** (Mie / TGF): multilayer, complex permittivity
  (absorption), magnetic layers, PEC limit, diffraction and surface-source
  ("antenna") problems, linear and circular polarization.
- **One solver per geometry**, all behind a uniform API:

  ===============  =====================================  ==========================
  Geometry         Solver                                 Regime
  ===============  =====================================  ==========================
  Sphere           ``SphereSolver``                       exact (Mie / TGF)
  Spheroid         ``SpheroidSolver``                     quasi-static (Rayleigh)
  Ellipsoid        ``EllipsoidSolver``                    quasi-static (Rayleigh)
  Cylinder (∞)     ``CylinderSolver``                     exact 2-D
  Cone             ``ConeSolver`` → sphere cluster + GMM  rigorous via decomposition
  Cluster          ``Cluster``                            self-consistent GMM
  ===============  =====================================  ==========================

- **Unified T-matrix interface** (spherical VSWF) + **GMM** assembly engine with
  closed-form Cruzan translation-addition coefficients (no dense-packing limit).
- **Honest scope.** Full-wave branches that require special-function machinery not
  available in SciPy (confocal spheroid, finite cylinder, sphero-conal cone) raise
  ``NotImplementedError`` rather than returning unverified results; complex bodies
  are handled rigorously by decomposition into the analytic sphere family + GMM.

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

   # Cluster of scatterers (GMM)
   s1 = gt.SphereSolver(0.5, [2.25], position=(-1, 0, 0)).as_scatterer()
   s2 = gt.SphereSolver(0.5, [2.25], position=(+1, 0, 0)).as_scatterer()
   gt.Cluster([s1, s2]).cross_sections(k=3.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=6)

   # Complex geometry: cone -> non-overlapping spheres -> GMM
   cone = gt.ConeSolver(apex=[0, 0, 0], axis=[0, 0, 1],
                        half_angle=0.5, height=8.0, eps=[2.25])
   scatterers, centers, radius = cone.decompose(spacing=2.0)
   gt.Cluster(scatterers).cross_sections(k=1.0, khat=(1, 0, 0), pol=(0, 0, 1), nmax=3)

See ``examples/`` for end-to-end scripts (Luneburg/Maxwell lenses, complex geometry).

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
     mie_core.py        importable mirror of the sphere core
     tmatrix.py         diagonal T-matrix + cross sections
     vswf.py            VSWF, Wigner-3j/Gaunt, closed-form translation (Cruzan)
     scatterer.py       Scatterer protocol + LayeredSphere
     gmm.py             Generalized Multiparticle Mie assembly engine
     ellipsoid.py,
     spheroid.py,
     cylinder.py,
     cone.py            per-geometry analytic solvers
     decompose.py       pack a body into non-overlapping spheres
     solvers.py         unified public API (SphereSolver, ... , Cluster, solve_*)
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

MIT — see `LICENSE <LICENSE>`_.

Citation
--------

If you use GreenTensor in academic work, please cite the IEEE papers above.
