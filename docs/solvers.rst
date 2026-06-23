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
   * - Cylinder (infinite)
     - :class:`~green_tensor.solvers.CylinderSolver`
     - :func:`~green_tensor.solvers.solve_cylinder`
     - **Exact** 2-D
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

:class:`~green_tensor.solvers.CylinderSolver` solves the infinite circular
cylinder at normal incidence exactly (TM/TE 2-D Mie series). The finite cylinder
is non-separable; ``finite()`` raises ``NotImplementedError`` — model finite bodies
with :class:`~green_tensor.solvers.ConeSolver`/``decompose`` or a prolate spheroid.

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
