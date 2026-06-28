# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""T-матрица в сферическом базисе ВСВФ — общий интерфейс всех примитивов.

Каждый рассеиватель (слоистая сфера; кластер сфер) описывается своей T-матрицей в
ЕДИНОМ базисе сферических векторных волновых функций, связывающей коэффициенты
падающего и рассеянного полей: ``c = T·d`` (см. GreenTensor_Theory.tex, разделы
«Постановка» и «Суперпозиция T-матриц»). Это объект, который потребляет движок сборки
кластера сфер (GMM, ``gmm.py``).

Конвенция Мищенко/Уотермана, время ``e^{-iωt}``, исходящая волна ~ ``h_n^{(1)}``.
Для СФЕРЫ T-матрица диагональна и не зависит от азимута ``m``:

    t^{(M, магнитн./TE)}_n = -b_n ,   t^{(N, электр./TM)}_n = -a_n ,

где ``a_n, b_n`` — стандартные коэффициенты Ми. Эффективности (нормировка на πR²,
параметр размера ``x = k0·R``):

    Q_sca  = (2/x²) Σ (2n+1) (|t^N_n|² + |t^M_n|²)
    Q_ext  = -(2/x²) Σ (2n+1) Re(t^M_n + t^N_n)
    Q_back = (1/x²) |Σ (2n+1) (-1)^n (t^M_n - t^N_n)|²

Внутри сечения считаются через a_n = -t^N_n, b_n = -t^M_n (формулы Bohren & Huffman),
чтобы исключить ошибки знака.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class DiagonalTMatrix:
    """Диагональная (по n, не зависящая от m) T-матрица — случай сферы.

    n   : порядки мультиполей 1..nmax;
    t_M : магнитные (TE) диагональные элементы, t^M_n = -b_n;
    t_N : электрические (TM) диагональные элементы, t^N_n = -a_n;
    x   : параметр размера k0·R для нормировки сечений на πR².
    """

    n: np.ndarray
    t_M: np.ndarray
    t_N: np.ndarray
    x: float

    def __post_init__(self):
        self.n = np.asarray(self.n)
        self.t_M = np.asarray(self.t_M, dtype=complex)
        self.t_N = np.asarray(self.t_N, dtype=complex)
        if not (len(self.n) == len(self.t_M) == len(self.t_N)):
            raise ValueError("n, t_M, t_N must have equal length")
        if not (self.x > 0 and np.isfinite(self.x)):
            raise ValueError("x must be a positive finite size parameter")

    # коэффициенты Ми в стандартной конвенции B&H
    @property
    def a(self) -> np.ndarray:  # электрический (N)
        return -self.t_N

    @property
    def b(self) -> np.ndarray:  # магнитный (M)
        return -self.t_M

    # ------------------------------------------------------------- сечения
    def q_sca(self) -> float:
        w = 2 * self.n + 1
        return float((2.0 / self.x**2) * np.sum(w * (np.abs(self.a) ** 2 + np.abs(self.b) ** 2)))

    def q_ext(self) -> float:
        w = 2 * self.n + 1
        return float((2.0 / self.x**2) * np.sum(w * np.real(self.a + self.b)))

    def q_abs(self) -> float:
        return self.q_ext() - self.q_sca()

    def q_back(self) -> float:
        w = 2 * self.n + 1
        s = np.sum(w * ((-1) ** self.n) * (self.a - self.b))
        return float((1.0 / self.x**2) * np.abs(s) ** 2)

    def cross_sections(self) -> dict:
        return {"q_sca": self.q_sca(), "q_ext": self.q_ext(),
                "q_abs": self.q_abs(), "q_back": self.q_back()}


def from_ab(n, a, b, x: float) -> DiagonalTMatrix:
    """Построить диагональную T-матрицу из коэффициентов Ми a_n (электр.), b_n (магн.)."""
    return DiagonalTMatrix(n=np.asarray(n), t_M=-np.asarray(b, dtype=complex),
                           t_N=-np.asarray(a, dtype=complex), x=float(x))


def from_mie_coeffs(Mn, Nn, x: float) -> DiagonalTMatrix:
    """Построить T-матрицу из коэффициентов канонического ядра сферы (01_sphere/sphere_core).

    ВНИМАНИЕ по именованию (footgun для GMM): в ядре ``Mn`` получается из импеданса Z и
    физически является ЭЛЕКТРИЧЕСКИМ коэффициентом a_n (TM), а ``Nn`` — из адмитанса Y и
    является МАГНИТНЫМ b_n (TE) (см. теорию §«Слоистая сфера», eq:Mn-Nn). Поэтому
    электрический блок t_N = −Mn (= −a_n), магнитный блок t_M = −Nn (= −b_n).
    """
    Mn = np.asarray(Mn, dtype=complex)
    Nn = np.asarray(Nn, dtype=complex)
    n = np.arange(1, len(Mn) + 1)
    return DiagonalTMatrix(n=n, t_M=-Nn, t_N=-Mn, x=float(x))   # t_M=−b_n(магн), t_N=−a_n(электр)


def sphere_tmatrix(mie) -> DiagonalTMatrix:
    """Адаптер: каноническая сфера (sphere_core.MieSphere) -> диагональная T-матрица.

    Принимает любой объект с методом ``coefficients() -> (Mn, Nn)`` и атрибутом ``k0``.
    """
    Mn, Nn = mie.coefficients()
    return from_mie_coeffs(Mn, Nn, mie.k0)
