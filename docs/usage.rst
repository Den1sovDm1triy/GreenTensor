.. SPDX-License-Identifier: MIT
.. Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

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

Ellipsoid and spheroid (quasi-static)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # triaxial ellipsoid, random-orientation average
   gt.solve_ellipsoid(a=2, b=1, c=1, eps=3 + 0.1j, k=0.5)

   # field along a chosen principal axis
   el = gt.EllipsoidSolver(2, 1, 1, eps=3 + 0.1j)
   el.cross_sections(k=0.5, axis=0)
   el.depolarization_factors()         # (L_a, L_b, L_c), sum = 1

   # spheroid in the Rayleigh limit = ellipsoid with b = a
   gt.EllipsoidSolver(1, 1, 2, eps=3 + 0.1j).cross_sections(k=0.5)

Rigorous non-spherical primitives (EBCM/TGF)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # full-wave spheroid / finite cylinder / cone — assembled through Cluster (GMM)
   sph = gt.Spheroid(position=(0, 0, 0), a_eq=0.4, c_ax=0.8, eps=2.25)
   gt.Cluster([sph]).cross_sections(k=2.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=7)

   cyl = gt.FiniteCylinder(position=(0, 0, 0), radius=0.35, half_length=0.5, eps=2.25)
   cone = gt.Cone(position=(0, 0, 0), radius=0.35, height=0.8, eps=2.25,
                  euler=(0.0, 0.5, 0.0))     # arbitrary orientation (Wigner-D)
   gt.Cluster([cyl, cone]).cross_sections(k=2.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=7)

   # layered primitive (Peterson–Ström recursion): eps inner->outer + a_norm boundaries
   gt.Spheroid((0, 0, 0), 0.4, 0.8, eps=[4.0, 2.25], a_norm=[0.6, 1.0])

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

   # cone -> non-overlapping spheres -> GMM
   cone = gt.ConeSolver(apex=[0, 0, 0], axis=[0, 0, 1],
                        half_angle=0.5, height=8.0, eps=[2.25])
   scatterers, centers, radius = cone.decompose(spacing=2.0)
   gt.Cluster(scatterers).cross_sections(k=1.0, khat=(1, 0, 0), pol=(0, 0, 1), nmax=3)

   # finite cylinder -> non-overlapping spheres -> GMM
   cyl = gt.FiniteCylinderSolver(center=[0, 0, 0], radius=2.0, half_length=4.0,
                                 eps=[2.25], axis=2)
   scatterers, centers, radius = cyl.decompose(spacing=1.0)
   gt.Cluster(scatterers).cross_sections(k=0.8, khat=(1, 0, 0), pol=(0, 0, 1), nmax=3)

Not-implemented (honest) branches
---------------------------------

The full-wave spheroid, finite cylinder and cone are now provided rigorously by the
EBCM primitives above. The remaining honest stub is the legacy 2-D module's finite
branch, which points to the rigorous :class:`~green_tensor.scatterer.FiniteCylinder`:

.. code-block:: python

   gt.CylinderSolver(1, 2.0).finite()   # raises NotImplementedError -> use gt.FiniteCylinder

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
