# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""scatterer — общий контракт рассеивателя для движка сборки GMM.

Каждый примитив реализует ``Scatterer``: положение, описывающий радиус и
T-матрицу в сферическом базисе ВСВФ (вектор для GMM). Сфера — каноническая
реализация поверх ``01_sphere.py``; остальные геометрии
(сфероид/эллипсоид/цилиндр/конус) добавляются как новые классы с тем же контрактом.

Вектор T для GMM хранит диагональ (M,N): c^M = t_M·a^M, c^N = t_N·a^N.
Для сферы t_M = −b_n (магн.), t_N = −a_n (электр.) — см. tmatrix.py.
"""
from __future__ import annotations

import math
import warnings
from typing import Protocol, runtime_checkable

import numpy as np

from . import ebcm
from . import sphere_core
from . import vswf


@runtime_checkable
class Scatterer(Protocol):
    position: np.ndarray

    def bounding_radius(self) -> float: ...

    def t_vector(self, k: float, nmax: int) -> np.ndarray:
        """Диагональный вектор T длины 2K (блоки M, затем N) для мод mode_list(nmax).

        Несферические примитивы вместо этого реализуют :meth:`t_matrix` (полная
        2K×2K матрица); gmm._full_t выбирает t_matrix, если он есть.
        """
        ...


def _is_layered(eps) -> bool:
    """eps задаёт несколько слоёв (список/массив) или один однородный материал (скаляр)?"""
    return np.iterable(eps) and not isinstance(eps, (str, bytes))


class LayeredSphere:
    """Радиально-слоистая сфера (каноническое ядро) как рассеиватель GMM."""

    def __init__(self, position, radius: float, eps, a_norm=None, miy=None):
        self.position = np.asarray(position, dtype=float)
        self.radius = float(radius)
        self.eps = list(eps)
        # нормированные радиусы границ слоёв (внешний = 1); по умолчанию однородная
        self.a_norm = list(a_norm) if a_norm is not None else [1.0] * len(self.eps)
        self.miy = list(miy) if miy is not None else None

    def bounding_radius(self) -> float:
        return self.radius

    def t_vector(self, k: float, nmax: int) -> np.ndarray:
        x = k * self.radius
        toch = max(nmax + 2, int(math.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0)))
        mie = sphere_core.MieSphere(k0=x, a=self.a_norm, eps=self.eps,
                                    miy=self.miy, toch=toch)
        Mn, Nn = mie.coefficients()
        modes = vswf.mode_list(nmax)
        tM = np.array([-Nn[n - 1] for (n, m) in modes], dtype=complex)  # −b_n
        tN = np.array([-Mn[n - 1] for (n, m) in modes], dtype=complex)  # −a_n
        return np.concatenate([tM, tN])


# --------------------------------------------------------------------------- #
# Несферические EBCM-примитивы как рассеиватели GMM (полная T-матрица)
# Тело строится в собственной системе (ось симметрии ∥ z); произвольное ПОЛОЖЕНИЕ —
# трансляция GMM; произвольная ОРИЕНТАЦИЯ — углы Эйлера euler=(α,β,γ), z-y-z, активный
# поворот R=Rz(α)Ry(β)Rz(γ): T_global = D(α,β,γ)·T·D^† (vswf.rotate_tmatrix). Для
# осесимметричного тела γ (поворот вокруг своей оси) — пустая операция.
# --------------------------------------------------------------------------- #
def _oriented(T, nmax, euler):
    """Применить ориентацию (углы Эйлера) к T-матрице тела (если поворот не нулевой)."""
    a, b, g = euler
    if a == 0.0 and b == 0.0 and g == 0.0:
        return T
    return vswf.rotate_tmatrix(T, nmax, a, b, g)


class Spheroid:
    """Сфероид как рассеиватель GMM через EBCM. eps скаляр — однородный; eps список
    (изнутри наружу) + a_norm — слоистый (mu аналогично). euler=(α,β,γ) — ориентация
    оси симметрии (по умолчанию ∥ z)."""

    def __init__(self, position, a_eq: float, c_ax: float, eps, mu=1.0, *,
                 a_norm=None, euler=(0.0, 0.0, 0.0)):
        self.position = np.asarray(position, dtype=float)
        self.a_eq = float(a_eq)
        self.c_ax = float(c_ax)
        self.eps = eps
        self.mu = mu
        self.a_norm = a_norm
        self.euler = tuple(float(x) for x in euler)

    def bounding_radius(self) -> float:
        return max(self.a_eq, self.c_ax)

    def t_matrix(self, k: float, nmax: int) -> np.ndarray:
        curve = ebcm.spheroid_curve(self.a_eq, self.c_ax)
        if _is_layered(self.eps):
            mu = self.mu if _is_layered(self.mu) else [self.mu] * len(self.eps)
            if self.a_norm is None:
                raise ValueError("слоистый сфероид требует a_norm (границы слоёв, внешний=1)")
            T = ebcm.tmatrix_axisym_layered(curve, list(self.eps), list(mu),
                                            list(self.a_norm), k, nmax)
        else:
            T = ebcm.tmatrix_axisym(curve, k, complex(self.eps), complex(self.mu), nmax)
        return _oriented(T, nmax, self.euler)


class FiniteCylinder:
    """Конечный круговой цилиндр (ось ∥ z, радиус R, полудлина H) — EBCM-рассеиватель GMM.

    Острые рёбра 90° (сингулярность Мейкснера) разложение по сферическим ВСВФ представляет
    плохо: сходимость по nmax крайне медленная, T-матрица теряет ВЗАИМНОСТЬ (~10–30% ошибки),
    причём расширенная точность и квадратура НЕ помогают (это усечение, не обусловленность).
    Поэтому по умолчанию рёбра СКРУГЛЯЮТСЯ суперэллипсом степени ``edge_p`` (гладкая образующая
    ebcm.rounded_cylinder_curve) — EBCM сходится быстро и взаимность ~1e-2. ``edge_p=None`` —
    прежний острый цилиндр (кромочно-ограниченный, для совместимости). Слоистый путь
    (a_norm) — через острые сегменты (скругление слоёв пока не поддержано)."""

    def __init__(self, position, radius: float, half_length: float, eps, mu=1.0, *,
                 a_norm=None, euler=(0.0, 0.0, 0.0), edge_p: float | None = 6.0):
        self.position = np.asarray(position, dtype=float)
        self.radius = float(radius)
        self.half_length = float(half_length)
        self.eps = eps
        self.mu = mu
        self.a_norm = a_norm
        self.euler = tuple(float(x) for x in euler)
        self.edge_p = None if edge_p is None else float(edge_p)

    def bounding_radius(self) -> float:
        return math.hypot(self.radius, self.half_length)

    def t_matrix(self, k: float, nmax: int) -> np.ndarray:
        R, H = self.radius, self.half_length
        if _is_layered(self.eps):
            mu = self.mu if _is_layered(self.mu) else [self.mu] * len(self.eps)
            if self.a_norm is None:
                raise ValueError("слоистый цилиндр требует a_norm")
            # ЧЕСТНЫЙ GUARD: слоистая рекурсия Петерсона–Стрёма на ОСТРЫХ сегментах цилиндра
            # численно неустойчива (исходящие h_n на внутр. границах + рёбра 90° ⇒ нарушение
            # унитарности в разы). Скругление слоёв пока не реализовано. Для слоистого цилиндра
            # используйте бесконечный ТФГ-решатель LayeredCylinderSolver или однослойный
            # (скруглённый) FiniteCylinder.
            warnings.warn(
                "Слоистый конечный цилиндр (EBCM) численно НЕУСТОЙЧИВ (рёбра + рекурсия слоёв): "
                "T-матрица может грубо нарушать унитарность. Используйте LayeredCylinderSolver "
                "(бесконечный ТФГ) или однослойный FiniteCylinder.",
                stacklevel=2,
            )
            n_per = 2 * nmax + 6
            nphi = 2 * nmax + 2

            def builder(a):
                return ebcm.surface_segments(ebcm.cylinder_segments(a * R, a * H), n_per, nphi)
            T = ebcm.tmatrix_layered(builder, list(self.eps), list(mu),
                                     list(self.a_norm), k, nmax)
        elif self.edge_p is None:
            T = ebcm.tmatrix_axisym_segments(ebcm.cylinder_segments(R, H),
                                             k, complex(self.eps), complex(self.mu), nmax)
        else:
            T = ebcm.tmatrix_axisym(ebcm.rounded_cylinder_curve(R, H, p=self.edge_p),
                                    k, complex(self.eps), complex(self.mu), nmax)
        return _oriented(T, nmax, self.euler)


class Cone:
    """Конечный круговой конус (ось ∥ z, базовый радиус R, высота L) — EBCM-рассеиватель GMM.
    Острая вершина + ребро основания; для умеренных конусов EBCM сходится хорошо."""

    def __init__(self, position, radius: float, height: float, eps, mu=1.0, *,
                 a_norm=None, euler=(0.0, 0.0, 0.0)):
        self.position = np.asarray(position, dtype=float)
        self.radius = float(radius)
        self.height = float(height)
        self.eps = eps
        self.mu = mu
        self.a_norm = a_norm
        self.euler = tuple(float(x) for x in euler)

    def bounding_radius(self) -> float:
        # вершина z_a=3L/4, дно z=−L/4, радиус R: дальняя точка — вершина или ребро основания
        return max(0.75 * self.height, math.hypot(self.radius, 0.25 * self.height))

    def t_matrix(self, k: float, nmax: int) -> np.ndarray:
        R, L = self.radius, self.height
        if _is_layered(self.eps):
            mu = self.mu if _is_layered(self.mu) else [self.mu] * len(self.eps)
            if self.a_norm is None:
                raise ValueError("слоистый конус требует a_norm")
            n_per = 2 * nmax + 6
            nphi = 2 * nmax + 2

            def builder(a):
                return ebcm.surface_segments(ebcm.cone_segments(a * R, a * L), n_per, nphi)
            T = ebcm.tmatrix_layered(builder, list(self.eps), list(mu),
                                     list(self.a_norm), k, nmax)
        else:
            T = ebcm.tmatrix_axisym_segments(ebcm.cone_segments(R, L),
                                             k, complex(self.eps), complex(self.mu), nmax)
        return _oriented(T, nmax, self.euler)
