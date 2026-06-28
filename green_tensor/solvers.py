# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""solvers — единый публичный API: по одному решателю на геометрию.
solvers — unified public API: one solver per geometry.

RU: Каждый решатель экспонируется как самостоятельный класс с единообразным
интерфейсом (``cross_sections`` и т.п.) поверх уже проверенных ядер пакета. Классы
НЕ содержат собственной математики — только связывают проверенные функции
(01_sphere.py/tmatrix/cylinder/gmm) в удобный объект. Для краткого вызова есть
функции-обёртки ``solve_*``.

EN: Each solver is exposed as a standalone class with a uniform interface
(``cross_sections`` etc.) on top of the package's already-verified cores. The classes
carry NO math of their own — they only wire the verified functions
(01_sphere.py/tmatrix/cylinder/gmm) into a convenient object. Thin ``solve_*`` wrapper
functions are provided for one-line calls.

ОБЛАСТЬ ОХВАТА: только ТОЧНЫЕ аналитические решения метода ТФГ — геометрии с
разделением переменных (сфера, бесконечный цилиндр) и их строгая суперпозиция
(кластер невзаимопересекающихся сфер, GMM, теорема сложения Крузана–Стейна).
Приближённые/численные методы (EBCM для тел с рёбрами, квазистатика Рэлея,
разложение в кластер сфер) намеренно НЕ входят в библиотеку.
SCOPE: only EXACT analytic TGF solutions — separable geometries (sphere, infinite
cylinder) and their rigorous superposition (cluster of non-overlapping spheres, GMM,
Cruzan–Stein addition theorem). Approximate/numerical methods (EBCM for edged bodies,
Rayleigh quasi-statics, sphere-cluster decomposition) are intentionally OUT of scope.

Карта решателей / solver map (see GreenTensor_Theory.tex):
  * :class:`SphereSolver`    — слоистая сфера, точное ядро (Mie/TFG) / layered sphere, exact core;
  * :class:`CylinderSolver`  — бесконечный цилиндр, точная 2D-аналитика / infinite cylinder, exact 2D;
  * :class:`LayeredCylinderSolver` — слоистый/косой бесконечный цилиндр (ТФГ) / layered/oblique infinite cylinder (TGF);
  * :class:`Cluster`         — сборка кластера сфер (GMM) / assembly of a sphere cluster (GMM).

Конвенция / convention: e^{-iωt}, внешняя среда — вакуум / vacuum host, исходящая
волна / outgoing wave ~ h_n^{(1)}. ``k`` — волновое число во внешней среде / wavenumber
in the host; параметр размера / size parameter x = k·(радиус / radius).
"""
from __future__ import annotations

import numpy as np

from . import cylinder as _cylinder
from . import gmm as _gmm
from . import scatterer as _scatterer
from . import sphere_core as _sphere_core
from . import tmatrix as _tmatrix

__all__ = [
    "SphereSolver",
    "CylinderSolver", "LayeredCylinderSolver", "Cluster",
    "solve_sphere",
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
        """RU: Конечный цилиндр несепарабелен — НЕ реализован; точная аналитика существует
        только для БЕСКОНЕЧНОГО цилиндра (этот класс / :class:`LayeredCylinderSolver`).
        EN: The finite cylinder is non-separable — NOT implemented; exact analytics exist only
        for the INFINITE cylinder (this class / :class:`LayeredCylinderSolver`)."""
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
# Кластер — движок сборки (GMM) для набора сфер
# Cluster — assembly engine (GMM) for a set of spheres
# --------------------------------------------------------------------------- #
class Cluster:
    """RU: Самосогласованная сборка кластера сфер (GMM, теорема сложения Крузана–Стейна).
    Принимает объекты протокола :class:`green_tensor.scatterer.Scatterer`
    (``SphereSolver(...).as_scatterer()``). Описывающие сферы тел НЕ должны пересекаться.
    EN: Self-consistent assembly of a sphere cluster (GMM, Cruzan–Stein addition theorem).
    Accepts objects implementing the :class:`green_tensor.scatterer.Scatterer` protocol
    (``SphereSolver(...).as_scatterer()``). Circumscribing spheres must not overlap.
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
