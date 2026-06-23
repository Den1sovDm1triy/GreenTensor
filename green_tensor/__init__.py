"""GreenTensor — analytic electromagnetic scattering on a layered-sphere core.

RU: Библиотека аналитического анализа рассеяния электромагнитных волн. Математическим
ядром служит точное решение для радиально-слоистой сферы (метод тензорных функций
Грина / Ми, ``01_sphere.py`` / :mod:`green_tensor.mie_core`). Поверх ядра построено
семейство аналитических решателей (сфероид, эллипсоид, цилиндр, конус) с единым
T-матричным интерфейсом в сферическом базисе ВСВФ и единый движок сборки сложной
геометрии (GMM, :mod:`green_tensor.gmm`).

EN: A library for analytic electromagnetic-scattering analysis. Its mathematical core
is the exact solution for a radially layered sphere (tensor Green's function / Mie
method, ``01_sphere.py`` / :mod:`green_tensor.mie_core`). On top of the core a family
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

__version__ = "0.2.0"

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
    spheroid,
    tmatrix,
    vswf,
)

# Публичный API / public API
from .scatterer import LayeredSphere, Scatterer  # noqa: F401
from .tmatrix import DiagonalTMatrix  # noqa: F401
from .solvers import (  # noqa: F401
    Cluster,
    ConeSolver,
    CylinderSolver,
    EllipsoidSolver,
    SphereSolver,
    SpheroidSolver,
    solve_cluster,
    solve_cylinder,
    solve_ellipsoid,
    solve_sphere,
    solve_spheroid,
)

__all__ = [
    "__version__",
    # protocols / data types
    "Scatterer", "LayeredSphere", "DiagonalTMatrix",
    # solver classes
    "SphereSolver", "EllipsoidSolver", "SpheroidSolver",
    "CylinderSolver", "ConeSolver", "Cluster",
    # functional facade
    "solve_sphere", "solve_ellipsoid", "solve_spheroid",
    "solve_cylinder", "solve_cluster",
    # modules
    "mie_core", "tmatrix", "scatterer", "vswf", "gmm",
    "ellipsoid", "spheroid", "cylinder", "cone", "decompose", "solvers",
]
