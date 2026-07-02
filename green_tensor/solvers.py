# SPDX-License-Identifier: MIT

"""solvers — единый публичный API: по одному решателю на геометрию.

Каждый решатель экспонируется как самостоятельный класс с единообразным
интерфейсом (``cross_sections`` и т.п.) поверх уже проверенных ядер пакета. Классы
НЕ содержат собственной математики — только связывают проверенные функции
(01_sphere.py/tmatrix/cylinder/gmm) в удобный объект. Для краткого вызова есть
функции-обёртки ``solve_*``.

ОБЛАСТЬ ОХВАТА: только ТОЧНЫЕ аналитические решения метода ТФГ — геометрии с
разделением переменных (сфера, бесконечный цилиндр) и их строгая суперпозиция
(кластер невзаимопересекающихся сфер, GMM, теорема сложения Крузана–Стейна).
Приближённые/численные методы (EBCM для тел с рёбрами, квазистатика Рэлея,
разложение в кластер сфер) намеренно НЕ входят в библиотеку.

Карта решателей (подробности — GreenTensor_Theory.tex):
  * :class:`SphereSolver`    — слоистая сфера, точное ядро (Ми/ТФГ);
  * :class:`CylinderSolver`  — бесконечный однородный цилиндр, точная 2D-аналитика;
  * :class:`LayeredCylinderSolver` — слоистый/косой бесконечный цилиндр (ТФГ);
  * :class:`Cluster`         — сборка кластера непересекающихся сфер (GMM).

Конвенция: e^{-iωt}, внешняя среда — вакуум, исходящая волна ~ h_n^{(1)};
``k`` — волновое число во внешней среде, параметр размера x = k·радиус.
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
# --------------------------------------------------------------------------- #
class SphereSolver:
    """Радиально-слоистая сфера — точный решатель поверх ``01_sphere.py``.

    radius : внешний радиус;
    eps    : комплексные проницаемости слоёв изнутри наружу;
    a_norm : нормированные радиусы границ слоёв (внешний = 1); по умолчанию однородная;
    miy    : магнитные проницаемости слоёв (по умолчанию 1);
    position : центр сферы для сборки кластера."""

    def __init__(self, radius, eps, *, a_norm=None, miy=None, position=(0.0, 0.0, 0.0)):
        self.radius = float(radius)
        self.eps = list(eps)
        self.a_norm = list(a_norm) if a_norm is not None else [1.0] * len(self.eps)
        self.miy = list(miy) if miy is not None else None
        self.position = np.asarray(position, dtype=float)

    def mie(self, k: float, *, toch: int | None = None, **mie_kwargs) -> "_sphere_core.MieSphere":
        """Объект канонического ядра ``01_sphere.py`` для x = k·radius."""
        return _sphere_core.MieSphere(k0=k * self.radius, a=self.a_norm, eps=self.eps,
                                      miy=self.miy, toch=toch, **mie_kwargs)

    def cross_sections(self, k: float, *, toch: int | None = None) -> dict:
        """Эффективности Q_sca, Q_ext, Q_abs, Q_back (норм. на πR²)."""
        return self.mie(k, toch=toch).cross_sections()

    def t_matrix(self, k: float, *, toch: int | None = None) -> "_tmatrix.DiagonalTMatrix":
        """Диагональная T-матрица сферы в сферическом базисе ВСВФ."""
        return _tmatrix.sphere_tmatrix(self.mie(k, toch=toch))

    def pattern(self, k: float, theta=None, **mie_kwargs) -> dict:
        """Диаграмма рассеяния (линейная/круговая поляризация) — формулы ``01_sphere.py``."""
        return self.mie(k, **mie_kwargs).pattern(theta)

    def as_scatterer(self) -> "_scatterer.LayeredSphere":
        """Представление для движка сборки кластера (GMM)."""
        return _scatterer.LayeredSphere(self.position, self.radius, self.eps,
                                        a_norm=self.a_norm, miy=self.miy)


# --------------------------------------------------------------------------- #
# Бесконечный круговой цилиндр — точная 2D-аналитика
# --------------------------------------------------------------------------- #
class CylinderSolver:
    """Бесконечный круговой цилиндр (нормальное падение), точное разделение переменных.

    radius — радиус; eps — проницаемость тела; eps_m — проницаемость среды.
    Относительный показатель преломления m = √(eps/eps_m)."""

    def __init__(self, radius, eps, *, eps_m: complex = 1.0):
        self.radius = float(radius)
        self.eps = complex(eps)
        self.eps_m = complex(eps_m)

    @property
    def m(self) -> complex:
        """Относительный комплексный показатель преломления √(eps/eps_m)."""
        return np.sqrt(self.eps / self.eps_m)

    def cross_sections(self, k: float, *, mode: str = "TM", nmax: int | None = None) -> dict:
        """Эффективности Q_sca, Q_ext, Q_abs (норм. на диаметр). mode: 'TM' | 'TE'."""
        return _cylinder.cross_sections_infinite(self.m, k * self.radius, mode=mode, nmax=nmax)

    def finite(self, *args, **kwargs):
        """Конечный цилиндр несепарабелен — НЕ реализован; точная аналитика существует
        только для БЕСКОНЕЧНОГО цилиндра (этот класс / :class:`LayeredCylinderSolver`)."""
        return _cylinder.finite(*args, **kwargs)


class LayeredCylinderSolver:
    """Радиально-СЛОИСТЫЙ бесконечный цилиндр (нормальное падение), метод ТФГ
    (эквивалентные линии передачи, Дайлис–Шабунин 2017). Точная многослойная
    рекурсия; при одном слое совпадает с :class:`CylinderSolver`.

    radius : внешний радиус;
    eps    : проницаемости слоёв изнутри наружу;
    mu     : магнитные проницаемости слоёв (по умолчанию 1);
    a_norm : нормированные радиусы границ слоёв (внешний = 1, по возрастанию);
             по умолчанию — оболочки равной толщины. Внешняя среда — вакуум."""

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
        """Эффективности Q_sca, Q_ext, Q_abs (норм. на диаметр). theta — угол падения
        от оси (по умолчанию нормальное); при нормальном — mode 'TM'(E∥)/'TE'(H∥), при
        косом — TM-падение (вертикальная поляризация) со связью поляризаций."""
        x = k * self.radius
        if theta is None or abs(theta - np.pi / 2) < 1e-12:
            return _cylinder.cross_sections_layered(self.eps, self.mu, self.a_norm,
                                                   x, mode=mode, nmax=nmax)
        return _cylinder.cross_sections_oblique(self.eps, self.mu, self.a_norm,
                                               x, theta, nmax=nmax)

    def coeff(self, k: float, n: int, *, mode: str = "TM") -> complex:
        """Коэффициент рассеяния a_n гармоники n (нормальное падение)."""
        return _cylinder.layered_coeff(self.eps, self.mu, self.a_norm,
                                      k * self.radius, n, mode=mode)

    def coeff_oblique(self, k: float, theta: float, n: int):
        """Рассеянные (A_E со-пол, A_H кросс-пол) при косом TM-падении, гармоника n."""
        return _cylinder.layered_coeff_oblique(self.eps, self.mu, self.a_norm,
                                              k * self.radius, theta, n)


# --------------------------------------------------------------------------- #
# Кластер — движок сборки (GMM) для набора сфер
# --------------------------------------------------------------------------- #
class Cluster:
    """Самосогласованная сборка кластера сфер (GMM, теорема сложения Крузана–Стейна).
    Принимает объекты протокола :class:`green_tensor.scatterer.Scatterer`
    (``SphereSolver(...).as_scatterer()``). Описывающие сферы тел НЕ должны пересекаться."""

    def __init__(self, scatterers):
        self.scatterers = list(scatterers)

    def solve(self, k: float, khat, pol, nmax: int):
        """Решить кластер: (a, c, d) — возбуждающие/рассеянные/падающие коэффициенты."""
        return _gmm.solve_cluster(self.scatterers, k, khat, pol, nmax)

    def cross_sections(self, k: float, khat, pol, nmax: int) -> dict:
        """Сечения кластера C_sca (дальнее поле), C_ext (опт. теорема), C_abs."""
        _, c, _ = _gmm.solve_cluster(self.scatterers, k, khat, pol, nmax)
        return _gmm.cross_sections(self.scatterers, c, k, khat, pol)

    def scattered_field(self, c, k: float, pts) -> np.ndarray:
        """Полное рассеянное поле кластера в точках pts по коэффициентам c (из :meth:`solve`)."""
        return _gmm.scattered_field(self.scatterers, c, k, pts)


# --------------------------------------------------------------------------- #
# Тонкий функциональный фасад — короткий вызов «решить и получить сечения»
# --------------------------------------------------------------------------- #
def solve_sphere(radius, eps, k: float, *, a_norm=None, miy=None, toch: int | None = None) -> dict:
    """Сечения слоистой сферы одной строкой (см. :class:`SphereSolver`)."""
    return SphereSolver(radius, eps, a_norm=a_norm, miy=miy).cross_sections(k, toch=toch)


def solve_cylinder(radius, eps, k: float, *, eps_m: complex = 1.0,
                   mode: str = "TM", nmax: int | None = None) -> dict:
    """Эффективности бесконечного цилиндра (см. :class:`CylinderSolver`)."""
    return CylinderSolver(radius, eps, eps_m=eps_m).cross_sections(k, mode=mode, nmax=nmax)


def solve_layered_cylinder(radius, eps, k: float, *, mu=None, a_norm=None,
                           mode: str = "TM", theta: float | None = None,
                           nmax: int | None = None) -> dict:
    """Эффективности слоистого цилиндра (ТФГ); theta — угол падения от оси
    (косое падение со связью поляризаций); см. :class:`LayeredCylinderSolver`."""
    return LayeredCylinderSolver(radius, eps, mu=mu, a_norm=a_norm).cross_sections(
        k, mode=mode, theta=theta, nmax=nmax)


def solve_cluster(scatterers, k: float, khat, pol, nmax: int) -> dict:
    """Сечения кластера непересекающихся сфер (см. :class:`Cluster`)."""
    return Cluster(scatterers).cross_sections(k, khat, pol, nmax)
