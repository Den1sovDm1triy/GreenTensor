.. SPDX-License-Identifier: MIT

Usage
=====

.. _installation:

Installation
------------

.. code-block:: console

   $ pip install .            # from a checkout
   $ pip install ".[docs]"    # + Sphinx docs
   $ pip install ".[dev]"     # + pytest, coverage

Requires Python ≥ 3.9 with NumPy, SciPy, Matplotlib. The distribution is named
``greentensor``; the import package is ``green_tensor``.

Conventions
-----------

- Time dependence :math:`e^{-i\omega t}`; outgoing waves :math:`\sim h_n^{(1)}`.
- Vacuum host medium (:math:`\varepsilon_m = 1` unless ``eps_m`` is given).
- ``k`` is the wavenumber in the host; the size parameter is :math:`x = k\,R`.
- ``Im(eps) > 0`` denotes absorption; ``a_norm`` are layer-boundary radii
  normalized so the outer radius is 1.

Quick start
-----------

Layered sphere (exact core)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import green_tensor as gt

   # one-line facade: cross sections normalized to πR²
   gt.solve_sphere(radius=1.0, eps=[2.25], k=3.0)
   # -> {'q_sca': ..., 'q_ext': ..., 'q_abs': ..., 'q_back': ...}

   # object API
   sphere = gt.SphereSolver(radius=1.0, eps=[2.25])
   sphere.cross_sections(k=3.0)
   T = sphere.t_matrix(k=3.0)          # DiagonalTMatrix in the VSWF basis
   diagram = sphere.pattern(k=3.0)     # scattering diagram (linear/circular)

   # multilayer + absorption + magnetic layers
   gt.SphereSolver(radius=1.0, eps=[4 + 0.5j, 2.25],
                   a_norm=[0.5, 1.0], miy=[1.0, 1.0]).cross_sections(k=2.0)

Infinite cylinder (exact 2-D)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # homogeneous
   gt.solve_cylinder(radius=1.0, eps=4.0, k=2.0, mode="TM")   # or mode="TE"

   # radially layered (TGF), normal incidence
   gt.solve_layered_cylinder(radius=1.0, eps=[4 + 0.5j, 2.0], mu=[1, 1],
                             a_norm=[0.6, 1.0], k=2.0, mode="TM")

   # oblique incidence (theta measured from the cylinder axis; couples polarizations)
   import numpy as np
   gt.solve_layered_cylinder(radius=1.0, eps=[4 + 0.5j, 2.0], a_norm=[0.6, 1.0],
                             k=2.0, theta=np.deg2rad(45))
   # per-harmonic co/cross-polarized coefficients:
   gt.LayeredCylinderSolver(1.0, [4 + 0.5j, 2.0], a_norm=[0.6, 1.0]) \
     .coeff_oblique(k=2.0, theta=np.deg2rad(45), n=1)   # -> (A_E, A_H)

Clusters and complex geometry (GMM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # two coupled spheres
   s1 = gt.SphereSolver(0.5, [2.25], position=(-1, 0, 0)).as_scatterer()
   s2 = gt.SphereSolver(0.5, [2.25], position=(+1, 0, 0)).as_scatterer()
   cluster = gt.Cluster([s1, s2])
   cluster.cross_sections(k=3.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=6)
   a, c, d = cluster.solve(k=3.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=6)

   # a multi-sphere chain (circumscribing spheres must not overlap)
   spheres = [gt.SphereSolver(0.5, [2.25], position=(0, 0, z)).as_scatterer()
              for z in (-1.5, 0.0, 1.5)]
   gt.Cluster(spheres).cross_sections(k=3.0, khat=(1, 0, 0), pol=(0, 0, 1), nmax=6)

Out-of-scope (honest) branches
------------------------------

Non-separable geometries with no exact closed-form solution are out of scope and
raise :class:`NotImplementedError` rather than returning an unverified stub — e.g.
the finite cylinder:

.. code-block:: python

   gt.CylinderSolver(1, 2.0).finite()   # raises NotImplementedError (use the infinite cylinder)

Running the tests
-----------------

.. code-block:: console

   $ python3 tests/run_all.py     # self-contained, no pytest required
   $ pytest tests/                # equivalent

Browser front-end
-----------------

A small web UI for the layered-sphere core lives in ``webapp/``:

.. code-block:: console

   $ python3 webapp/server.py
