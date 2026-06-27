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
from typing import Protocol, runtime_checkable

import numpy as np

from . import sphere_core
from . import vswf


@runtime_checkable
class Scatterer(Protocol):
    position: np.ndarray

    def bounding_radius(self) -> float: ...

    def t_vector(self, k: float, nmax: int) -> np.ndarray:
        """Диагональный вектор T длины 2K (блоки M, затем N) для мод mode_list(nmax)."""
        ...


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
