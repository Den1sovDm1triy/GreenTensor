.. SPDX-License-Identifier: MIT
.. Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

Solvers
=======

GreenTensor exposes **one solver per geometry**, all behind a uniform object API
in :mod:`green_tensor.solvers` (re-exported at the top level as ``green_tensor.*``).
Each class wraps an already-verified mathematical core; it carries no math of its
own. Thin ``solve_*`` functions provide one-line access.

.. list-table:: Solver map
   :header-rows: 1
   :widths: 18 30 26 26

   * - Geometry
     - Class
     - Facade
     - Regime
   * - Sphere (layered)
     - :class:`~green_tensor.solvers.SphereSolver`
     - :func:`~green_tensor.solvers.solve_sphere`
     - **Exact** (Mie / TGF)
   * - Spheroid
     - :class:`~green_tensor.solvers.SpheroidSolver`
     - :func:`~green_tensor.solvers.solve_spheroid`
     - Quasi-static (Rayleigh)
   * - Ellipsoid (triaxial)
     - :class:`~green_tensor.solvers.EllipsoidSolver`
     - :func:`~green_tensor.solvers.solve_ellipsoid`
     - Quasi-static (Rayleigh)
   * - Cylinder (infinite, homogeneous)
     - :class:`~green_tensor.solvers.CylinderSolver`
     - :func:`~green_tensor.solvers.solve_cylinder`
     - **Exact** 2-D
   * - Cylinder (infinite, layered, TGF)
     - :class:`~green_tensor.solvers.LayeredCylinderSolver`
     - :func:`~green_tensor.solvers.solve_layered_cylinder`
     - **Exact** (normal + oblique)
   * - Cylinder (finite)
     - :class:`~green_tensor.solvers.FiniteCylinderSolver`
     - via :class:`~green_tensor.solvers.Cluster`
     - Rigorous (sphere decomposition)
   * - Cone
     - :class:`~green_tensor.solvers.ConeSolver`
     - via :class:`~green_tensor.solvers.Cluster`
     - Rigorous (sphere decomposition)
   * - Arbitrary cluster
     - :class:`~green_tensor.solvers.Cluster`
     - :func:`~green_tensor.solvers.solve_cluster`
     - Self-consistent GMM

Sphere — the exact core
-----------------------

:class:`~green_tensor.solvers.SphereSolver` is the canonical, *exact* solver
(tensor Green's function / Mie recursion). It supports radial multilayers, complex
permittivity (absorption), magnetic layers, the PEC limit, and both the
``diffraction`` and surface-source (``antenna``) problems.

Key methods: ``cross_sections(k)`` (efficiencies normalized to :math:`\pi R^2`),
``t_matrix(k)`` (a :class:`~green_tensor.tmatrix.DiagonalTMatrix` in the VSWF
basis), ``pattern(k)`` (scattering diagram), and ``as_scatterer()`` (to drop the
sphere into a :class:`~green_tensor.solvers.Cluster`).

Ellipsoid and spheroid — quasi-static
-------------------------------------

:class:`~green_tensor.solvers.EllipsoidSolver` and
:class:`~green_tensor.solvers.SpheroidSolver` solve the Rayleigh
(:math:`k a \ll 1`) limit exactly: geometric depolarization factors, homogeneous
and confocal-coated polarizabilities, Rayleigh cross sections, and the
electric-dipole T-matrix element. ``cross_sections(k)`` returns the
random-orientation average by default, or the response along a chosen principal
axis with ``axis=i``. The spheroid additionally offers closed-form depolarization
factors (``depolarization_factors(closed=True)``).

The full-wave confocal spheroid (separation of variables in the spheroidal basis
with complex ``c``) is **not implemented** — ``SpheroidSolver.full_wave()`` raises
``NotImplementedError``; SciPy provides only real-``c`` spheroidal functions.

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

The **finite** cylinder is non-separable (no direct full-wave solver). Use
:class:`~green_tensor.solvers.FiniteCylinderSolver`, which ``decompose``\ s the body
into a non-overlapping sphere cluster solved by :class:`~green_tensor.solvers.Cluster`
(GMM) — the same rigorous-analytic route as :class:`~green_tensor.solvers.ConeSolver`.
A prolate spheroid is an alternative quasi-static approximation.

Finite cylinder — rigorous via decomposition
---------------------------------------------

:class:`~green_tensor.solvers.FiniteCylinderSolver` (centre, radius, half-length,
axis) has no direct full-wave solver; ``decompose(spacing)`` packs the cylinder with
non-overlapping spheres (verified non-overlapping and inside the body) which feed the
spherical-VSWF :class:`~green_tensor.solvers.Cluster` (GMM). ``full_wave()`` raises
``NotImplementedError``.

.. warning::

   **Metals are rejected.** A non-overlapping sphere packing leaves vacuum gaps; for a
   metal body these gaps are open cavities that produce spurious inter-sphere
   re-reflections and resonances absent in the solid body (a metallic "sponge", not a
   solid). ``decompose`` raises ``ValueError`` when the **outermost** layer is metallic
   — ``Re(eps) < 0`` or skin depth ``δ = 1/(k·Im√(εμ)) < r_sphere`` (pass ``k`` to
   enable the skin-depth test). A metal **core fully enclosed by a dielectric outer
   layer** is allowed (gaps touch only dielectric; the metal is a closed volume inside
   each sphere). Use ``allow_metal=True`` to deliberately model a cluster of metal
   spheres. The same guard applies to :class:`~green_tensor.solvers.ConeSolver`. For a
   solid metal body use a single closed body (:class:`~green_tensor.solvers.SphereSolver`
   with a metal layer) or a surface method (EBCM/MoM).

Decomposition optimizers
~~~~~~~~~~~~~~~~~~~~~~~~~~

The packing leaves vacuum gaps, so the cluster's effective permittivity is **lower**
than the solid body's (a systematic bias). ``decompose(..., effective_medium=True)``
applies a closed-form **Maxwell–Garnett** correction: it chooses the sphere permittivity
so the effective medium of (spheres at filling fraction *f* in vacuum) equals the body's
true ε. In the quasi-static limit this is exact — the corrected sub-spheres reproduce the
solid body's polarizability identically (``Σα_corrected == α_solid``), whereas the
uncorrected packing is biased by a factor *f*.

Helpers: :func:`~green_tensor.decompose.maxwell_garnett_eps` /
:func:`~green_tensor.decompose.maxwell_garnett_effective`,
:func:`~green_tensor.decompose.packing_report` (filling, overlap margin, GMM size).

.. note::

   The correction is **filling-limited**: the maximum reachable ε is
   ``(1 + 2f)/(1 − f)`` (an infinite-ε inclusion would be metal, which is forbidden). A
   sparse cubic packing gives ``f ≈ 0.1–0.4`` (reaches ε up to ≈ 1.4–3), so MG only
   corrects modest-ε bodies unless the packing is made denser; it raises ``ValueError``
   when *f* is too low for the requested ε. It is also quasi-static (valid for
   ``k·r_sphere ≪ 1``): it corrects the effective bulk permittivity, not surface detail.

   Pass ``lattice="fcc"`` to ``decompose`` (face-centered cubic) for a ~√2× denser,
   still non-overlapping packing — higher filling *f* and thus a higher reachable ε for
   the Maxwell–Garnett correction, plus better boundary conformity. ``lattice="cubic"``
   (default) is the simple cubic lattice.

Cone — rigorous via decomposition
---------------------------------

The direct sphero-conal solver needs non-integer-degree Legendre functions and is
**not implemented** (``ConeSolver.full_wave()`` raises ``NotImplementedError``).
The rigorous, *analytic* route is ``ConeSolver.decompose(spacing)`` → a cluster of
non-overlapping spheres solved by GMM. This pattern (``decompose`` →
:class:`~green_tensor.solvers.Cluster`) handles arbitrary "brick-like" geometry
from the strictly analytic sphere family.

Cluster — the GMM assembly engine
----------------------------------

:class:`~green_tensor.solvers.Cluster` composes any list of objects implementing
the :class:`~green_tensor.scatterer.Scatterer` protocol into a single
self-consistent solution via the translation-addition theorems (closed-form
Cruzan coefficients — valid at any separation). ``cross_sections(...)`` returns the
far-field ``C_sca``, the optical-theorem ``C_ext`` and ``C_abs``; ``solve(...)``
returns the raw exciting/scattered/incident coefficients; ``scattered_field(c, k,
pts)`` evaluates the total scattered field.

See :doc:`api` for the complete signatures.
