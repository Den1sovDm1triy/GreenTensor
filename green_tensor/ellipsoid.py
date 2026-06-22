"""ellipsoid — квазистатический (рэлеевский) решатель для трёхосного эллипсоида.

Полноволновая векторная задача для эллипсоида не разделяется (функции Ламе
нетабулируемы выше малых степеней) — см. GreenTensor_Theory.tex, раздел «Эллипсоид».
Практический аналитический решатель существует в квазистатическом пределе (k·a ≪ 1):
факторы деполяризации L_i, поляризуемость α_i (однородный и конфокальный покрытый
эллипсоид, Bohren & Huffman / Sihvola–Lindell) и рэлеевские сечения. Электрический
диполь даёт вклад n=1 в T-матрицу для сборки GMM.

Конвенция: внешняя среда вакуум (ε_m=1); полуоси a≥b≥c; α_i — вдоль главной оси i.
"""
from __future__ import annotations

import numpy as np
from scipy import integrate


def depolarization_factors(a: float, b: float, c: float):
    """Геометрические факторы деполяризации L_a, L_b, L_c (L_a+L_b+L_c=1)."""
    abc = a * b * c

    def Li(ai):
        f = lambda q: 1.0 / ((ai**2 + q) * np.sqrt((a**2 + q) * (b**2 + q) * (c**2 + q)))
        val, _ = integrate.quad(f, 0.0, np.inf)
        return 0.5 * abc * val

    return Li(a), Li(b), Li(c)


def polarizability_homogeneous(a, b, c, eps, eps_m: complex = 1.0):
    """Поляризуемость α_i однородного эллипсоида вдоль трёх главных осей."""
    V = (4.0 / 3.0) * np.pi * a * b * c
    L = depolarization_factors(a, b, c)
    de = eps - eps_m
    return np.array([V * de / (eps_m + Li * de) for Li in L], dtype=complex)


def polarizability_coated_confocal(a2, b2, c2, f: float, eps1, eps2, eps_m: complex = 1.0):
    """Поляризуемость конфокального ПОКРЫТОГО эллипсоида (Bohren & Huffman 5.34).

    a2,b2,c2 — полуоси внешней (мантии) поверхности; f=V_core/V_shell — доля объёма ядра;
    eps1 — ядро, eps2 — мантия. Конфокальное ядро имеет полуоси (a2²−ξ)^½ и т.д.
    """
    V2 = (4.0 / 3.0) * np.pi * a2 * b2 * c2
    L2 = depolarization_factors(a2, b2, c2)
    # полуоси конфокального ядра: a1²=a2²−d, b1²=b2²−d, c1²=c2²−d с d из доли объёма f
    # объём ядра = f·V2 ⇒ a1 b1 c1 = f·a2 b2 c2 ; находим d численно
    d = _confocal_offset(a2, b2, c2, f)
    a1, b1, c1 = np.sqrt(a2**2 - d), np.sqrt(b2**2 - d), np.sqrt(c2**2 - d)
    L1 = depolarization_factors(a1, b1, c1)
    alpha = np.empty(3, dtype=complex)
    for i in range(3):
        num = (eps2 - eps_m) * (eps2 + (eps1 - eps2) * (L1[i] - f * L2[i])) \
            + f * eps2 * (eps1 - eps2)
        den = (eps2 + (eps1 - eps2) * (L1[i] - f * L2[i])) * (eps_m + (eps2 - eps_m) * L2[i]) \
            + f * L2[i] * eps2 * (eps1 - eps2)
        alpha[i] = V2 * num / den
    return alpha


def _confocal_offset(a2, b2, c2, f: float) -> float:
    """Найти смещение d: (a2²−d)(b2²−d)(c2²−d) = f²·(a2 b2 c2)² (объём ядра=f·V)."""
    target = (f * a2 * b2 * c2) ** 2
    lo, hi = 0.0, c2**2 * (1 - 1e-12)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        val = (a2**2 - mid) * (b2**2 - mid) * (c2**2 - mid)
        if val > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def rayleigh_cross_sections(alpha_eff: complex, k: float) -> dict:
    """Рэлеевские сечения для эффективной поляризуемости α вдоль поля E."""
    c_sca = (k**4) / (6.0 * np.pi) * abs(alpha_eff) ** 2
    c_abs = k * np.imag(alpha_eff)
    return {"c_sca": c_sca, "c_abs": c_abs, "c_ext": c_sca + c_abs}


def dipole_t_scalar(alpha: complex, k: float) -> complex:
    """Электрич.-дипольный элемент T_N(n=1) из скалярной поляризуемости (рэлеевский предел)."""
    return 1j * k**3 * alpha / (6.0 * np.pi)
