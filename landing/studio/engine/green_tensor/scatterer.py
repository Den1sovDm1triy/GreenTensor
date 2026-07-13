# SPDX-License-Identifier: MIT

"""scatterer — общий контракт рассеивателя для движка сборки GMM.

Рассеиватель реализует ``Scatterer``: положение, описывающий радиус и T-матрицу в
сферическом базисе ВСВФ (диагональный вектор для GMM). Единственный примитив —
радиально-слоистая сфера :class:`LayeredSphere` (каноническое точное ядро
``01_sphere.py``); набор таких сфер связывается строгой теоремой сложения в
:mod:`green_tensor.gmm` (кластер невзаимопересекающихся сфер).

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
        """Диагональный вектор T длины 2K (блоки M, затем N) для мод mode_list(nmax).

        Рассеиватель с недиагональной T-матрицей может вместо этого реализовать
        :meth:`t_matrix` (полная 2K×2K матрица); gmm._full_t выбирает t_matrix, если он есть.
        """
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
        # Ядро сферы хранит СОПРЯЖЁННЫЕ коэффициенты Ми: Mn = a_n*, Nn = b_n*
        # (знак мнимой части инвертируется на последнем шаге расчёта, см. 01_sphere.py,
        # calculate_scattering_coefficients). T-матрица в конвенции пакета e^{-iωt}
        # строится из истинных a_n, b_n: t_M = −b_n = −Nn*, t_N = −a_n = −Mn*.
        # Без сопряжения фазы рассеянного поля неверны: нарушается граничное условие
        # на поверхности (тест PEC-сферы) и связь тел в кластере GMM; сечения
        # одиночной сферы к сопряжению нечувствительны, поэтому ловится только
        # полевыми и кластерными проверками.
        a_n = np.conj(Mn)
        b_n = np.conj(Nn)
        modes = vswf.mode_list(nmax)
        tM = np.array([-b_n[n - 1] for (n, m) in modes], dtype=complex)
        tN = np.array([-a_n[n - 1] for (n, m) in modes], dtype=complex)
        return np.concatenate([tM, tN])
