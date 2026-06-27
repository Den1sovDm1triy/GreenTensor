# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""solvers — единый публичный API: по одному решателю на геометрию.
solvers — unified public API: one solver per geometry.

RU: Каждый примитив семейства экспонируется как самостоятельный класс-решатель с
единообразным интерфейсом (``cross_sections`` и т.п.) поверх уже проверенных ядер
пакета. Классы НЕ содержат собственной математики — только связывают проверенные
функции (01_sphere.py/tmatrix/ellipsoid/cylinder/decompose/gmm) в удобный объект.
Для краткого вызова есть функции-обёртки ``solve_*``.

EN: Each primitive of the family is exposed as a standalone solver class with a
uniform interface (``cross_sections`` etc.) on top of the package's already-verified
cores. The classes carry NO math of their own — they only wire the verified functions
(01_sphere.py/tmatrix/ellipsoid/cylinder/decompose/gmm) into a convenient object.
Thin ``solve_*`` wrapper functions are provided for one-line calls.

Карта решателей / solver map (see GreenTensor_Theory.tex):
  * :class:`SphereSolver`    — слоистая сфера, точное ядро (Mie/TFG) / layered sphere, exact core;
  * :class:`EllipsoidSolver` — трёхосный эллипсоид, квазистатика / triaxial ellipsoid, quasi-static;
  * :class:`CylinderSolver`  — бесконечный цилиндр, точная 2D-аналитика / infinite cylinder, exact 2D;
  * :class:`LayeredCylinderSolver` — слоистый/косой бесконечный цилиндр (ТФГ) / layered/oblique infinite cylinder (TGF);
  * :class:`ConeSolver`, :class:`FiniteCylinderSolver` — decompose-fallback (кластер сфер + GMM)
    / decompose fallback (sphere cluster + GMM);
  * :class:`Cluster`         — сборка набора рассеивателей (GMM) / assembly of scatterers (GMM).

RU: СТРОГИЕ полноволновые примитивы — ``green_tensor.{Spheroid, FiniteCylinder, Cone}``
(метод EBCM/ТФГ); собираются через :class:`Cluster` (GMM) и дают точные сечения
одиночного тела или кластера. Это основной путь для несферических тел.
EN: The RIGOROUS full-wave primitives are ``green_tensor.{Spheroid, FiniteCylinder, Cone}``
(EBCM/TGF method), assembled via :class:`Cluster` (GMM) for exact single-body or
cluster cross sections — the main route for non-spherical bodies.

Конвенция / convention: e^{-iωt}, внешняя среда — вакуум / vacuum host, исходящая
волна / outgoing wave ~ h_n^{(1)}. ``k`` — волновое число во внешней среде / wavenumber
in the host; параметр размера / size parameter x = k·(радиус / radius).
"""
from __future__ import annotations

import numpy as np

from . import cylinder as _cylinder
from . import decompose as _decompose
from . import ellipsoid as _ellipsoid
from . import gmm as _gmm
from . import scatterer as _scatterer
from . import sphere_core as _sphere_core
from . import tmatrix as _tmatrix

__all__ = [
    "SphereSolver", "EllipsoidSolver",
    "CylinderSolver", "LayeredCylinderSolver", "FiniteCylinderSolver",
    "ConeSolver", "Cluster",
    "solve_sphere", "solve_ellipsoid",
    "solve_cylinder", "solve_layered_cylinder", "solve_cluster",
]


# --------------------------------------------------------------------------- #
# Сфера — каноническое точное ядро (метод Ми / ТФГ)
# Sphere — canonical exact core (Mie / TFG method)
# --------------------------------------------------------------------------- #
class SphereSolver:
    """RU: Радиально-слоистая сфера — точный решатель поверх ``01_sphere.py``.
    EN: Radially layered sphere — exact solver on top of ``01_sphere.py``.

    radius : внешний радиус / outer radius;
    eps    : комплексные проницаемости слоёв изнутри наружу / layer permittivities, inner→outer;
    a_norm : нормированные радиусы границ (внешний = 1); по умолчанию однородная
             / normalized layer-boundary radii (outer = 1); default homogeneous;
    miy    : магнитные проницаемости слоёв (по умолчанию 1) / layer permeabilities (default 1);
    position : центр сферы для сборки кластера / sphere centre for cluster assembly.
    """

    def __init__(self, radius, eps, *, a_norm=None, miy=None, position=(0.0, 0.0, 0.0)):
        self.radius = float(radius)
        self.eps = list(eps)
        self.a_norm = list(a_norm) if a_norm is not None else [1.0] * len(self.eps)
        self.miy = list(miy) if miy is not None else None
        self.position = np.asarray(position, dtype=float)

    def mie(self, k: float, *, toch: int | None = None, **mie_kwargs) -> "_sphere_core.MieSphere":
        """RU: Объект канонического ядра ``01_sphere.py`` для x = k·radius.
        EN: Canonical ``01_sphere.py`` core object for x = k·radius."""
        return _sphere_core.MieSphere(k0=k * self.radius, a=self.a_norm, eps=self.eps,
                                      miy=self.miy, toch=toch, **mie_kwargs)

    def cross_sections(self, k: float, *, toch: int | None = None) -> dict:
        """RU: Эффективности Q_sca, Q_ext, Q_abs, Q_back (норм. на πR²).
        EN: Efficiencies Q_sca, Q_ext, Q_abs, Q_back (normalized to πR²)."""
        return self.mie(k, toch=toch).cross_sections()

    def t_matrix(self, k: float, *, toch: int | None = None) -> "_tmatrix.DiagonalTMatrix":
        """RU: Диагональная T-матрица сферы в сферическом базисе ВСВФ.
        EN: Diagonal sphere T-matrix in the spherical VSWF basis."""
        return _tmatrix.sphere_tmatrix(self.mie(k, toch=toch))

    def pattern(self, k: float, theta=None, **mie_kwargs) -> dict:
        """RU: Диаграмма рассеяния (линейная/круговая поляризация) — формулы ``01_sphere.py``.
        EN: Scattering pattern (linear/circular polarization) — ``01_sphere.py`` formulas."""
        return self.mie(k, **mie_kwargs).pattern(theta)

    def as_scatterer(self) -> "_scatterer.LayeredSphere":
        """RU: Представление для движка сборки кластера (GMM).
        EN: Representation for the cluster-assembly engine (GMM)."""
        return _scatterer.LayeredSphere(self.position, self.radius, self.eps,
                                        a_norm=self.a_norm, miy=self.miy)


# --------------------------------------------------------------------------- #
# Трёхосный эллипсоид — квазистатика (рэлеевский предел)
# Triaxial ellipsoid — quasi-static (Rayleigh) limit
# --------------------------------------------------------------------------- #
class EllipsoidSolver:
    """RU: Однородный трёхосный эллипсоид в квазистатическом пределе (k·a ≪ 1).
    EN: Homogeneous triaxial ellipsoid in the quasi-static limit (k·a ≪ 1).

    Полуоси / semi-axes a ≥ b ≥ c; eps — тело / body permittivity; eps_m — среда / host.
    """

    def __init__(self, a, b, c, eps, *, eps_m: complex = 1.0, position=(0.0, 0.0, 0.0)):
        self.a, self.b, self.c = float(a), float(b), float(c)
        self.eps = complex(eps)
        self.eps_m = complex(eps_m)
        self.position = np.asarray(position, dtype=float)

    def depolarization_factors(self):
        """RU: Факторы деполяризации (L_a, L_b, L_c), Σ = 1.
        EN: Depolarization factors (L_a, L_b, L_c), Σ = 1."""
        return _ellipsoid.depolarization_factors(self.a, self.b, self.c)

    def polarizability(self) -> np.ndarray:
        """RU: Поляризуемость α вдоль трёх главных осей.
        EN: Polarizability α along the three principal axes."""
        return _ellipsoid.polarizability_homogeneous(self.a, self.b, self.c,
                                                     self.eps, self.eps_m)

    def cross_sections(self, k: float, *, orientation_average: bool = True,
                       axis: int | None = None) -> dict:
        """RU: Рэлеевские сечения C_sca, C_abs, C_ext. orientation_average=True — среднее
        по случайной ориентации; axis=i — поле вдоль главной оси i (0/1/2).
        EN: Rayleigh cross sections C_sca, C_abs, C_ext. orientation_average=True —
        random-orientation average; axis=i — field along principal axis i (0/1/2)."""
        if axis is None and orientation_average:
            return _ellipsoid.orientation_average_cross_sections(
                self.a, self.b, self.c, self.eps, k, self.eps_m)
        i = 0 if axis is None else int(axis)
        return _ellipsoid.rayleigh_cross_sections(self.polarizability()[i], k)

    def dipole_t(self, k: float, *, axis: int = 0) -> complex:
        """RU: Электрич.-дипольный элемент T_N(n=1) вдоль оси axis (рэлеевский предел).
        EN: Electric-dipole element T_N(n=1) along ``axis`` (Rayleigh limit)."""
        return _ellipsoid.dipole_t_scalar(self.polarizability()[axis], k)


# --------------------------------------------------------------------------- #
# Сфероид (вытянутый / сплюснутый) — строгий полноволновой EBCM (ТФГ)
# Spheroid (prolate / oblate) — rigorous full-wave EBCM (TGF)
#   → публичный примитив green_tensor.Spheroid в green_tensor.Cluster (GMM):
#       gt.Cluster([gt.Spheroid(pos, a_eq, c_ax, eps)]).cross_sections(k, khat, pol, nmax)
#     Рэлеевский предел — через EllipsoidSolver(a_eq, a_eq, c_ax, …) (сфероид = эллипсоид b=a).
#   → public primitive green_tensor.Spheroid inside green_tensor.Cluster (GMM);
#     Rayleigh limit via EllipsoidSolver(a_eq, a_eq, c_ax, …) (spheroid = ellipsoid b=a).
# --------------------------------------------------------------------------- #


# --------------------------------------------------------------------------- #
# Бесконечный круговой цилиндр — точная 2D-аналитика
# Infinite circular cylinder — exact 2D analytics
# --------------------------------------------------------------------------- #
class CylinderSolver:
    """RU: Бесконечный круговой цилиндр (нормальное падение), точное разделение переменных.
    EN: Infinite circular cylinder (normal incidence), exact separation of variables.

    radius — радиус / radius; eps — проницаемость / permittivity; eps_m — среда / host.
    Относительный показатель / relative index m = √(eps/eps_m).
    """

    def __init__(self, radius, eps, *, eps_m: complex = 1.0):
        self.radius = float(radius)
        self.eps = complex(eps)
        self.eps_m = complex(eps_m)

    @property
    def m(self) -> complex:
        """RU: Относительный комплексный показатель преломления √(eps/eps_m).
        EN: Relative complex refractive index √(eps/eps_m)."""
        return np.sqrt(self.eps / self.eps_m)

    def cross_sections(self, k: float, *, mode: str = "TM", nmax: int | None = None) -> dict:
        """RU: Эффективности Q_sca, Q_ext, Q_abs (норм. на диаметр). mode: 'TM' | 'TE'.
        EN: Efficiencies Q_sca, Q_ext, Q_abs (normalized to diameter). mode: 'TM' | 'TE'."""
        return _cylinder.cross_sections_infinite(self.m, k * self.radius, mode=mode, nmax=nmax)

    def finite(self, *args, **kwargs):
        """RU: Конечный цилиндр (мод-матчинг/EBCM) — НЕ реализован (см. ConeSolver/decompose).
        EN: Finite cylinder (mode-matching/EBCM) — NOT implemented (see ConeSolver/decompose)."""
        return _cylinder.finite(*args, **kwargs)


class LayeredCylinderSolver:
    """RU: Радиально-СЛОИСТЫЙ бесконечный цилиндр (нормальное падение), метод ТФГ
    (эквивалентные линии передачи, Дайлис–Шабунин 2017). Точная многослойная
    рекурсия; при одном слое совпадает с :class:`CylinderSolver`.
    EN: Radially LAYERED infinite cylinder (normal incidence), TGF method
    (equivalent transmission lines, Daylis–Shabunin 2017). Exact multilayer
    recursion; reduces to :class:`CylinderSolver` for a single layer.

    radius : внешний радиус / outer radius;
    eps    : проницаемости слоёв изнутри наружу / layer permittivities, inner→outer;
    mu     : магнитные проницаемости слоёв (по умолчанию 1) / layer permeabilities (default 1);
    a_norm : нормированные радиусы границ (внешний = 1, по возрастанию); по умолчанию
             равные толщины / normalized boundary radii (outer = 1, increasing);
             default equal-thickness shells. Внешняя среда — вакуум / vacuum host.
    """

    def __init__(self, radius, eps, *, mu=None, a_norm=None):
        self.radius = float(radius)
        self.eps = list(eps)
        n = len(self.eps)
        self.mu = list(mu) if mu is not None else [1.0] * n
        if a_norm is not None:
            self.a_norm = list(a_norm)
        else:
            self.a_norm = [(i + 1) / n for i in range(n)]   # равные толщины, внешний = 1

    def cross_sections(self, k: float, *, mode: str = "TM", theta: float | None = None,
                       nmax: int | None = None) -> dict:
        """RU: Эффективности Q_sca, Q_ext, Q_abs (норм. на диаметр). theta — угол падения
        от оси (по умолчанию нормальное); при нормальном — mode 'TM'(E∥)/'TE'(H∥), при
        косом — TM-падение (вертикальная поляризация) со связью поляризаций.
        EN: Efficiencies Q_sca, Q_ext, Q_abs (per diameter). theta = incidence angle from
        the axis (default normal); at normal — mode 'TM'(E∥)/'TE'(H∥), at oblique —
        TM incidence (vertical polarization) with polarization coupling."""
        x = k * self.radius
        if theta is None or abs(theta - np.pi / 2) < 1e-12:
            return _cylinder.cross_sections_layered(self.eps, self.mu, self.a_norm,
                                                   x, mode=mode, nmax=nmax)
        return _cylinder.cross_sections_oblique(self.eps, self.mu, self.a_norm,
                                               x, theta, nmax=nmax)

    def coeff(self, k: float, n: int, *, mode: str = "TM") -> complex:
        """RU: Коэффициент рассеяния a_n гармоники n (нормальное падение).
        EN: scattering coefficient a_n of harmonic n (normal incidence)."""
        return _cylinder.layered_coeff(self.eps, self.mu, self.a_norm,
                                      k * self.radius, n, mode=mode)

    def coeff_oblique(self, k: float, theta: float, n: int):
        """RU: Рассеянные (A_E со-пол, A_H кросс-пол) при косом TM-падении, гармоника n.
        EN: scattered (A_E co-pol, A_H cross-pol) for oblique TM incidence, harmonic n."""
        return _cylinder.layered_coeff_oblique(self.eps, self.mu, self.a_norm,
                                              k * self.radius, theta, n)


# --------------------------------------------------------------------------- #
# Конус — decompose-fallback (кластер сфер + GMM). Строгий полноволновой конус —
# EBCM-примитив green_tensor.Cone в Cluster (см. ниже).
# Cone — decompose fallback (sphere cluster + GMM). Rigorous full-wave cone is the
# EBCM primitive green_tensor.Cone inside Cluster.
# --------------------------------------------------------------------------- #
class ConeSolver:
    """RU: Конечный круговой конус: вершина apex, ось axis, полуугол half_angle, высота height.
    Это decompose-fallback (разложение в кластер сфер → :class:`Cluster` (GMM)) для
    конусоподобных тел. СТРОГИЙ полноволновой конус — примитив ``green_tensor.Cone``
    (EBCM) в :class:`Cluster`.
    EN: Finite circular cone: apex, axis, half-angle, height. This is the decompose
    fallback (sphere-cluster → :class:`Cluster` (GMM)) for cone-like bodies. The
    RIGOROUS full-wave cone is the ``green_tensor.Cone`` primitive (EBCM) in :class:`Cluster`.
    """

    def __init__(self, apex, axis, half_angle: float, height: float, eps, *,
                 a_norm=None, miy=None):
        self.apex = np.asarray(apex, dtype=float)
        self.axis = np.asarray(axis, dtype=float)
        self.half_angle = float(half_angle)
        self.height = float(height)
        self.eps = list(eps)
        self.a_norm = a_norm
        self.miy = miy

    def indicator(self):
        """RU: Логический индикатор «точка внутри конуса».
        EN: Boolean indicator "point inside the cone"."""
        return _decompose.cone_indicator(self.apex, self.axis, self.half_angle, self.height)

    def decompose(self, spacing: float, *, fill: float = 0.45,
                  k: float | None = None, allow_metal: bool = False,
                  effective_medium: bool = False, lattice: str = "cubic"):
        """RU: Разложить в непересекающиеся сферы: (scatterers, centers, radius) для GMM.
        Металл с металлической внешней средой запрещён (k — для скин-слойной проверки,
        allow_metal=True — обход). effective_medium=True — поправка Максвелла–Гарнетта на ε
        (однородный диэлектрик, квазистатика). lattice='cubic'|'fcc' (ГЦК плотнее).
        EN: Decompose into non-overlapping spheres for GMM. Metal outer layer rejected
        (pass k for the skin-depth test; allow_metal=True to override). effective_medium=True
        applies the Maxwell–Garnett eps correction. lattice='cubic'|'fcc' (FCC is denser)."""
        return _decompose.decompose_cone(self.apex, self.axis, self.half_angle, self.height,
                                         spacing, self.eps, fill=fill,
                                         a_norm=self.a_norm, miy=self.miy,
                                         k=k, allow_metal=allow_metal,
                                         effective_medium=effective_medium, lattice=lattice)


# --------------------------------------------------------------------------- #
# Конечный цилиндр — строгая аналитика через разложение в кластер сфер + GMM
# Finite cylinder — rigorous analytics via sphere-cluster decomposition + GMM
# --------------------------------------------------------------------------- #
class FiniteCylinderSolver:
    """RU: Конечный круговой цилиндр (центр, радиус, полудлина, ось). Прямого
    полноволнового решателя нет (несепарабелен — см. :class:`CylinderSolver`.finite);
    строго аналитический путь — :meth:`decompose` в кластер непересекающихся сфер
    → :class:`Cluster` (GMM). Для БЕСКОНЕЧНОГО цилиндра используйте
    :class:`CylinderSolver` / :class:`LayeredCylinderSolver` (точная аналитика).
    EN: Finite circular cylinder (centre, radius, half-length, axis). No direct
    full-wave solver (non-separable — see :class:`CylinderSolver`.finite); the
    rigorous-analytic route is :meth:`decompose` into a non-overlapping sphere
    cluster → :class:`Cluster` (GMM). For an INFINITE cylinder use
    :class:`CylinderSolver` / :class:`LayeredCylinderSolver` (exact analytics).
    """

    def __init__(self, center, radius: float, half_length: float, eps, *,
                 axis: int = 2, a_norm=None, miy=None):
        self.center = np.asarray(center, dtype=float)
        self.radius = float(radius)
        self.half_length = float(half_length)
        self.eps = list(eps)
        self.axis = int(axis)
        self.a_norm = a_norm
        self.miy = miy

    def indicator(self):
        """RU: Логический индикатор «точка внутри цилиндра».
        EN: Boolean indicator "point inside the cylinder"."""
        return _decompose.cylinder_indicator(self.center, self.radius,
                                             self.half_length, self.axis)

    def decompose(self, spacing: float, *, fill: float = 0.45,
                  k: float | None = None, allow_metal: bool = False,
                  effective_medium: bool = False, lattice: str = "cubic"):
        """RU: Разложить в непересекающиеся сферы: (scatterers, centers, radius) для GMM.
        Металл с металлической внешней средой запрещён (k — для скин-слойной проверки,
        allow_metal=True — обход). effective_medium=True — поправка Максвелла–Гарнетта на ε
        (однородный диэлектрик, квазистатика). lattice='cubic'|'fcc' (ГЦК плотнее).
        EN: Decompose into non-overlapping spheres for GMM. Metal outer layer rejected
        (pass k for the skin-depth test; allow_metal=True to override). effective_medium=True
        applies the Maxwell–Garnett eps correction. lattice='cubic'|'fcc' (FCC is denser)."""
        return _decompose.decompose_cylinder(self.center, self.radius, self.half_length,
                                            spacing, self.eps, axis=self.axis, fill=fill,
                                            a_norm=self.a_norm, miy=self.miy,
                                            k=k, allow_metal=allow_metal,
                                            effective_medium=effective_medium, lattice=lattice)

    def full_wave(self, *args, **kwargs):
        """RU: Прямой решатель конечного цилиндра — НЕ реализован (мод-матчинг/EBCM).
        EN: Direct finite-cylinder solver — NOT implemented (mode-matching/EBCM)."""
        return _cylinder.finite(*args, **kwargs)


# --------------------------------------------------------------------------- #
# Кластер — движок сборки (GMM) для произвольного набора рассеивателей
# Cluster — assembly engine (GMM) for an arbitrary set of scatterers
# --------------------------------------------------------------------------- #
class Cluster:
    """RU: Самосогласованная сборка списка рассеивателей (GMM, трансляционные теоремы).
    Принимает любые объекты протокола :class:`green_tensor.scatterer.Scatterer`
    (например ``SphereSolver(...).as_scatterer()`` или результат ``ConeSolver.decompose``).
    EN: Self-consistent assembly of a list of scatterers (GMM, translation theorems).
    Accepts any object implementing the :class:`green_tensor.scatterer.Scatterer` protocol
    (e.g. ``SphereSolver(...).as_scatterer()`` or the result of ``ConeSolver.decompose``).
    """

    def __init__(self, scatterers):
        self.scatterers = list(scatterers)

    def solve(self, k: float, khat, pol, nmax: int):
        """RU: Решить кластер: (a, c, d) — возбуждающие/рассеянные/падающие коэффициенты.
        EN: Solve the cluster: (a, c, d) — exciting/scattered/incident coefficients."""
        return _gmm.solve_cluster(self.scatterers, k, khat, pol, nmax)

    def cross_sections(self, k: float, khat, pol, nmax: int) -> dict:
        """RU: Сечения кластера C_sca (дальнее поле), C_ext (опт. теорема), C_abs.
        EN: Cluster cross sections C_sca (far field), C_ext (optical theorem), C_abs."""
        _, c, _ = _gmm.solve_cluster(self.scatterers, k, khat, pol, nmax)
        return _gmm.cross_sections(self.scatterers, c, k, khat, pol)

    def scattered_field(self, c, k: float, pts) -> np.ndarray:
        """RU: Полное рассеянное поле кластера в точках pts по коэффициентам c (из :meth:`solve`).
        EN: Total scattered field of the cluster at ``pts`` from coefficients c (from :meth:`solve`)."""
        return _gmm.scattered_field(self.scatterers, c, k, pts)


# --------------------------------------------------------------------------- #
# Тонкий функциональный фасад — короткий вызов «решить и получить сечения»
# Thin functional facade — one-line "solve and get cross sections"
# --------------------------------------------------------------------------- #
def solve_sphere(radius, eps, k: float, *, a_norm=None, miy=None, toch: int | None = None) -> dict:
    """RU: Сечения слоистой сферы одной строкой / EN: layered-sphere cross sections in one line
    (см. / see :class:`SphereSolver`)."""
    return SphereSolver(radius, eps, a_norm=a_norm, miy=miy).cross_sections(k, toch=toch)


def solve_ellipsoid(a, b, c, eps, k: float, *, eps_m: complex = 1.0,
                    orientation_average: bool = True, axis: int | None = None) -> dict:
    """RU: Рэлеевские сечения трёхосного эллипсоида / EN: Rayleigh cross sections of a triaxial
    ellipsoid (см. / see :class:`EllipsoidSolver`)."""
    return EllipsoidSolver(a, b, c, eps, eps_m=eps_m).cross_sections(
        k, orientation_average=orientation_average, axis=axis)


def solve_cylinder(radius, eps, k: float, *, eps_m: complex = 1.0,
                   mode: str = "TM", nmax: int | None = None) -> dict:
    """RU: Эффективности бесконечного цилиндра / EN: infinite-cylinder efficiencies
    (см. / see :class:`CylinderSolver`)."""
    return CylinderSolver(radius, eps, eps_m=eps_m).cross_sections(k, mode=mode, nmax=nmax)


def solve_layered_cylinder(radius, eps, k: float, *, mu=None, a_norm=None,
                           mode: str = "TM", theta: float | None = None,
                           nmax: int | None = None) -> dict:
    """RU: Эффективности слоистого цилиндра (ТФГ); theta — угол падения от оси (косое
    падение со связью поляризаций) / EN: layered-cylinder efficiencies (TGF); theta =
    incidence angle from the axis (oblique incidence with polarization coupling)
    (см. / see :class:`LayeredCylinderSolver`)."""
    return LayeredCylinderSolver(radius, eps, mu=mu, a_norm=a_norm).cross_sections(
        k, mode=mode, theta=theta, nmax=nmax)


def solve_cluster(scatterers, k: float, khat, pol, nmax: int) -> dict:
    """RU: Сечения произвольного кластера рассеивателей / EN: cross sections of an arbitrary
    cluster of scatterers (см. / see :class:`Cluster`)."""
    return Cluster(scatterers).cross_sections(k, khat, pol, nmax)
