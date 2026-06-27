# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""GreenTensor — analytic electromagnetic scattering on a layered-sphere core.

RU: Библиотека аналитического анализа рассеяния электромагнитных волн. Математическим
ядром служит точное решение для радиально-слоистой сферы (метод тензорных функций
Грина / Ми, ``01_sphere.py``; импортируемый фасад — :mod:`green_tensor.sphere_core`).
Несферические тела (сфероид, конечный цилиндр, конус, в т.ч. слоистые) решаются СТРОГО
полноволновым методом EBCM/ТФГ — примитивы ``Spheroid``/``FiniteCylinder``/``Cone`` с
единым T-матричным интерфейсом в сферическом базисе ВСВФ; сложная геометрия собирается
единым движком GMM (:mod:`green_tensor.gmm`). Дополняют семейство замкнутые аналитики
(трёхосный эллипсоид — квазистатика; бесконечный/слоистый/косой цилиндр — точная 2D) и
fallback-разложение произвольных тел в кластер сфер (:mod:`green_tensor.decompose`).

EN: A library for analytic electromagnetic-scattering analysis. Its mathematical core
is the exact solution for a radially layered sphere (tensor Green's function / Mie
method, ``01_sphere.py``; importable facade :mod:`green_tensor.sphere_core`). Non-spherical
bodies (spheroid, finite cylinder, cone, incl. layered) are solved RIGOROUSLY by the
full-wave EBCM/TGF method — primitives ``Spheroid``/``FiniteCylinder``/``Cone`` sharing
one T-matrix interface in the spherical-VSWF basis; complex geometry is composed by a
single GMM engine (:mod:`green_tensor.gmm`). Complementary closed-form analytics (triaxial
ellipsoid — quasi-static; infinite/layered/oblique cylinder — exact 2D) and a fallback
decomposition of arbitrary bodies into a sphere cluster (:mod:`green_tensor.decompose`)
round out the family.

Быстрый старт / quick start::

    import green_tensor as gt

    # сфера / sphere
    gt.solve_sphere(radius=1.0, eps=[2.25], k=3.0)

    # строгий полноволновой несферический примитив (EBCM/ТФГ) / rigorous full-wave
    # non-spherical primitive (EBCM/TGF)
    cone = gt.Cone(position=(0, 0, 0), radius=0.35, height=0.8, eps=2.25)
    gt.Cluster([cone]).cross_sections(k=2.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=7)

    # кластер / cluster (GMM)
    s1 = gt.SphereSolver(0.5, [2.25], position=(-1, 0, 0)).as_scatterer()
    s2 = gt.SphereSolver(0.5, [2.25], position=(+1, 0, 0)).as_scatterer()
    gt.Cluster([s1, s2]).cross_sections(k=3.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=6)

См. / see GreenTensor_Theory.tex и docs/.
"""
from __future__ import annotations

__version__ = "0.4.0"

# Низкоуровневые модули / low-level modules
from . import (  # noqa: F401
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
    "ellipsoid", "cylinder", "decompose", "solvers",
]
