.. SPDX-License-Identifier: MIT
.. Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

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
