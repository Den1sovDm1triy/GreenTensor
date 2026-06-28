.. SPDX-License-Identifier: MIT

Solvers
=======

GreenTensor exposes only **exact analytic** tensor Green's function (TGF) solvers —
geometries that admit separation of variables, plus their rigorous superposition:

* the radially **layered sphere** (Mie / TGF closed form);
* the infinite **circular cylinder** (homogeneous and radially layered, normal and
  oblique incidence; exact 2-D analytics);
* a **cluster of non-overlapping spheres**, coupled self-consistently by the
  Generalized Multiparticle Mie (GMM) method via the Cruzan–Stein addition theorem.

Approximate or purely numerical methods (EBCM for edged bodies, Rayleigh
quasi-statics, decomposition of a body into packed spheres) are intentionally **not**
part of the library.

.. list-table:: Solver map
   :header-rows: 1
   :widths: 30 30 26 14

   * - Geometry
     - Class
     - Facade
     - Regime
   * - Sphere (layered)
     - :class:`~green_tensor.solvers.SphereSolver`
     - :func:`~green_tensor.solvers.solve_sphere`
     - **Exact** (Mie / TGF)
   * - Cylinder (infinite, homogeneous)
     - :class:`~green_tensor.solvers.CylinderSolver`
     - :func:`~green_tensor.solvers.solve_cylinder`
     - **Exact** 2-D
   * - Cylinder (infinite, layered, TGF)
     - :class:`~green_tensor.solvers.LayeredCylinderSolver`
     - :func:`~green_tensor.solvers.solve_layered_cylinder`
     - **Exact** (normal + oblique)
   * - Cluster of non-overlapping spheres
     - :class:`~green_tensor.solvers.Cluster`
     - :func:`~green_tensor.solvers.solve_cluster`
     - **Exact** superposition (GMM)

Sphere — the exact core
-----------------------

:class:`~green_tensor.solvers.SphereSolver` is the canonical, *exact* solver
(tensor Green's function / Mie recursion) built on the canonical kernel
``green_tensor/01_sphere.py`` (importable facade
:mod:`green_tensor.sphere_core`). It supports radial multilayers, complex
permittivity (absorption), magnetic layers and the PEC limit, in both the plane-wave
**diffraction** and surface-source (**antenna**) problems
(``pattern(k, problem="antenna")``).

Key methods: ``cross_sections(k)`` (efficiencies normalized to :math:`\pi R^2`),
``t_matrix(k)`` (a :class:`~green_tensor.tmatrix.DiagonalTMatrix` in the VSWF
basis), ``pattern(k)`` (scattering diagram), and ``as_scatterer()`` (to drop the
sphere into a :class:`~green_tensor.solvers.Cluster`).

Cylinder — exact 2-D
--------------------

:class:`~green_tensor.solvers.CylinderSolver` solves the infinite *homogeneous*
circular cylinder at normal incidence exactly (TM/TE 2-D Mie series).

:class:`~green_tensor.solvers.LayeredCylinderSolver` extends this to a radially
**multilayer** cylinder via the tensor Green's function (TGF) / equivalent
transmission-line method (Daylis–Shabunin), covering both **normal** and
**oblique** incidence. At normal incidence the TM/TE polarizations decouple into
independent scalar layer recursions; at oblique incidence (``theta`` ≠ π/2) they
couple (the ``n·k_z`` term), producing co- and cross-polarized scattering. It is
verified to machine precision against Bohren–Huffman (normal) and the Wait /
Kavaklıoğlu isolated-cylinder coefficients (oblique), reduces to the homogeneous
solver at one layer, and conserves energy.

The *finite* cylinder is non-separable and has no exact closed-form solution in this
basis, so it is out of scope; ``CylinderSolver.finite`` raises
:class:`NotImplementedError` rather than returning an unverified stub.

Cluster — the GMM assembly engine (non-overlapping spheres)
-----------------------------------------------------------

:class:`~green_tensor.solvers.Cluster` composes a list of spheres (objects
implementing the :class:`~green_tensor.scatterer.Scatterer` protocol, e.g.
``SphereSolver(...).as_scatterer()``) into a single self-consistent solution via the
translation-addition theorem (closed-form Cruzan coefficients — valid at any
separation). The **circumscribing spheres must not overlap** — that is the condition
under which the addition theorem (and hence the GMM superposition) is rigorous.

``cross_sections(...)`` returns the far-field ``C_sca``, the optical-theorem
``C_ext`` and ``C_abs``; ``solve(...)`` returns the raw
exciting/scattered/incident coefficients; ``scattered_field(c, k, pts)`` evaluates
the total scattered field.

See :doc:`api` for the complete signatures.
