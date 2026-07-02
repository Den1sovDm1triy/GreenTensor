# SPDX-License-Identifier: MIT

"""GreenTensor — analytic electromagnetic scattering on a layered-sphere core.

Библиотека аналитического анализа рассеяния электромагнитных волн. ОБЛАСТЬ —
только ТОЧНЫЕ аналитические решения метода тензорных функций Грина (ТФГ): геометрии
с разделением переменных и их строгая суперпозиция.

  * радиально-слоистая СФЕРА — точное матядро ``01_sphere.py`` (Ми/ТФГ; импортируемый
    фасад :mod:`green_tensor.sphere_core`);
  * бесконечный круговой ЦИЛИНДР — одно-/многослойный, нормальное и косое падение,
    точная 2D-аналитика (ТФГ, эквивалентные линии передачи; :mod:`green_tensor.cylinder`);
  * КЛАСТЕР невзаимопересекающихся сфер — строгая суперпозиция через теорему сложения
    Крузана–Стейна (:mod:`green_tensor.gmm`).

Приближённые и чисто численные методы (EBCM для тел с рёбрами, квазистатика Рэлея,
разложение тел в кластер сфер) в библиотеку намеренно НЕ входят — только верифицируемая
аналитика.

A library for analytic electromagnetic-scattering analysis. SCOPE — only EXACT
analytic tensor-Green's-function (TGF) solutions: separable geometries and their
rigorous superposition.

  * radially layered SPHERE — exact core ``01_sphere.py`` (Mie/TGF; importable facade
    :mod:`green_tensor.sphere_core`);
  * infinite circular CYLINDER — single/multilayer, normal and oblique incidence, exact
    2D analytics (TGF, equivalent transmission lines; :mod:`green_tensor.cylinder`);
  * CLUSTER of non-overlapping spheres — rigorous superposition via the Cruzan–Stein
    addition theorem (:mod:`green_tensor.gmm`).

Approximate and purely numerical methods (EBCM for edged bodies, Rayleigh quasi-statics,
decomposition of bodies into a sphere cluster) are intentionally NOT included — only
verifiable analytics.

Быстрый старт / quick start::

    import green_tensor as gt

    # сфера / sphere
    gt.solve_sphere(radius=1.0, eps=[2.25], k=3.0)

    # бесконечный (слоистый) цилиндр / infinite (layered) cylinder
    gt.solve_layered_cylinder(radius=1.0, eps=[2.25], k=2.0)

    # кластер сфер / sphere cluster (GMM)
    s1 = gt.SphereSolver(0.5, [2.25], position=(-1, 0, 0)).as_scatterer()
    s2 = gt.SphereSolver(0.5, [2.25], position=(+1, 0, 0)).as_scatterer()
    gt.Cluster([s1, s2]).cross_sections(k=3.0, khat=(0, 0, 1), pol=(1, 0, 0), nmax=6)

См. / see GreenTensor_Theory.tex и docs/.
"""
from __future__ import annotations

__version__ = "0.5.0"

# Низкоуровневые модули / low-level modules
from . import (  # noqa: F401
    cylinder,
    gmm,
    scatterer,
    solvers,
    sphere_core,
    tmatrix,
    vswf,
)

# Публичный API / public API
from .scatterer import (  # noqa: F401
    LayeredSphere,
    Scatterer,
)
from .tmatrix import DiagonalTMatrix  # noqa: F401
from .solvers import (  # noqa: F401
    Cluster,
    CylinderSolver,
    LayeredCylinderSolver,
    SphereSolver,
    solve_cluster,
    solve_cylinder,
    solve_layered_cylinder,
    solve_sphere,
)

__all__ = [
    "__version__",
    # protocols / data types
    "Scatterer", "LayeredSphere", "DiagonalTMatrix",
    # solver classes
    "SphereSolver",
    "CylinderSolver", "LayeredCylinderSolver", "Cluster",
    # functional facade
    "solve_sphere",
    "solve_cylinder", "solve_layered_cylinder", "solve_cluster",
    # modules
    "sphere_core", "tmatrix", "scatterer", "vswf", "gmm",
    "cylinder", "solvers",
]
