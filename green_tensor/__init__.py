# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""GreenTensor — analytic electromagnetic scattering on a layered-sphere core.

RU: Библиотека аналитического анализа рассеяния электромагнитных волн. Математическим
ядром служит точное решение для радиально-слоистой сферы (метод тензорных функций
Грина / Ми, ``01_sphere.py``; импортируемый фасад — :mod:`green_tensor.sphere_core`). Поверх ядра построено
семейство аналитических решателей (сфероид, эллипсоид, цилиндр, конус) с единым
T-матричным интерфейсом в сферическом базисе ВСВФ и единый движок сборки сложной
геометрии (GMM, :mod:`green_tensor.gmm`).

EN: A library for analytic electromagnetic-scattering analysis. Its mathematical core
is the exact solution for a radially layered sphere (tensor Green's function / Mie
method, ``01_sphere.py``; importable facade :mod:`green_tensor.sphere_core`). On top of the core a family
of analytic solvers (spheroid, ellipsoid, cylinder, cone) shares one T-matrix interface
in the spherical-VSWF basis, composed by a single complex-geometry assembly engine
(GMM, :mod:`green_tensor.gmm`).

Быстрый старт / quick start::

    import green_tensor as gt

    # сфера / sphere
    gt.solve_sphere(radius=1.0, eps=[2.25], k=3.0)

    # кластер / cluster (GMM)
    s1 = gt.SphereSolver(0.5, [2.25], position=(-1, 0, 0)).as_scatterer()
    s2 = gt.SphereSolver(0.5, [2.25], position=(+1, 0, 0)).as_scatterer()
    gt.Cluster([s1, s2]).cross_sections(k=3.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=6)

См. / see GreenTensor_Theory.tex и docs/.
"""
from __future__ import annotations

__version__ = "0.3.0"

# Низкоуровневые модули / low-level modules
from . import (  # noqa: F401
    cone,
    cylinder,
    decompose,
    ellipsoid,
    gmm,
    mie_core,
    scatterer,
    solvers,
    sphere_core,
    tmatrix,
    vswf,
)

# Публичный API / public API
# Строгие полноволновые EBCM-примитивы (ТФГ) — собираются через :class:`Cluster` (GMM).
# Rigorous full-wave EBCM primitives (TGF) — assembled via :class:`Cluster` (GMM).
from .scatterer import (  # noqa: F401
    Cone,
    FiniteCylinder,
    LayeredSphere,
    Scatterer,
    Spheroid,
)
from .tmatrix import DiagonalTMatrix  # noqa: F401
from .solvers import (  # noqa: F401
    Cluster,
    ConeSolver,
    CylinderSolver,
    EllipsoidSolver,
    FiniteCylinderSolver,
    LayeredCylinderSolver,
    SphereSolver,
    solve_cluster,
    solve_cylinder,
    solve_ellipsoid,
    solve_layered_cylinder,
    solve_sphere,
)

__all__ = [
    "__version__",
    # protocols / data types
    "Scatterer", "LayeredSphere", "DiagonalTMatrix",
    # строгие EBCM-примитивы (ТФГ) / rigorous EBCM primitives (TGF)
    "Spheroid", "FiniteCylinder", "Cone",
    # solver classes
    "SphereSolver", "EllipsoidSolver",
    "CylinderSolver", "LayeredCylinderSolver", "FiniteCylinderSolver",
    "ConeSolver", "Cluster",
    # functional facade
    "solve_sphere", "solve_ellipsoid",
    "solve_cylinder", "solve_layered_cylinder", "solve_cluster",
    # modules
    "mie_core", "sphere_core", "tmatrix", "scatterer", "vswf", "gmm",
    "ellipsoid", "cylinder", "cone", "decompose", "solvers",
]
