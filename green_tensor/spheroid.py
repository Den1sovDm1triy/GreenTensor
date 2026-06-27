# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""spheroid — вытянутый/сплюснутый сфероид.

КВАЗИСТАТИЧЕСКИЙ (рэлеевский) решатель реализован точно: сфероид — частный случай
эллипсоида с двумя равными полуосями, поэтому опирается на проверенную
деполяризационную машину (ellipsoid.py) и снабжён ЗАМКНУТЫМИ формулами факторов
деполяризации (Bohren & Huffman), сверяемыми с численным интегралом.

ПОЛНОВОЛНОВОЙ (конфокальный многослойный) сфероид — метод разделения переменных в
сфероидальном базисе (Asano–Yamamoto / Farafonov–Voshchinnikov, см.
GreenTensor_Theory.tex, раздел «Сфероид»). Он требует устойчивого вычисления
сфероидальных волновых функций с КОМПЛЕКСНЫМ параметром c=k·d/2 (поглощающие слои),
чего нет в scipy (pro_*/obl_* — только вещественный c). Эта ветка намеренно НЕ
реализована заглушкой-фейком: она поднимает NotImplementedError со ссылкой на
требуемую библиотеку, чтобы не выдавать непроверенный результат.

Геометрия: ось симметрии c_ax, экваториальная полуось a_eq. prolate: c_ax>a_eq,
oblate: c_ax<a_eq.
"""
from __future__ import annotations

import math

import numpy as np

from . import ellipsoid as _el


def depolarization(a_eq: float, c_ax: float):
    """Факторы деполяризации (L_eq, L_eq, L_ax) сфероида (численно, через эллипсоид)."""
    L_eq, _, L_ax = _el.depolarization_factors(a_eq, a_eq, c_ax)
    return L_eq, L_eq, L_ax


def depolarization_closed(a_eq: float, c_ax: float):
    """Замкнутая форма (Bohren & Huffman): возвращает (L_eq, L_ax)."""
    if abs(c_ax - a_eq) < 1e-12:
        return 1.0 / 3.0, 1.0 / 3.0
    if c_ax > a_eq:  # вытянутый (prolate), ось c — длинная
        e = math.sqrt(1.0 - (a_eq / c_ax) ** 2)
        L_ax = (1 - e**2) / e**2 * (-1.0 + (1.0 / (2 * e)) * math.log((1 + e) / (1 - e)))
    else:            # сплюснутый (oblate), ось c — короткая
        e = math.sqrt(1.0 - (c_ax / a_eq) ** 2)
        g = math.sqrt((1 - e**2) / e**2)
        L_ax = (1.0 / e**2) * (1.0 - g * math.asin(e))
    L_eq = 0.5 * (1.0 - L_ax)
    return L_eq, L_ax


def polarizability(a_eq: float, c_ax: float, eps, eps_m: complex = 1.0):
    """Поляризуемость (α_eq, α_eq, α_ax) однородного сфероида."""
    return _el.polarizability_homogeneous(a_eq, a_eq, c_ax, eps, eps_m)


def rayleigh_cross_sections(alpha_eff: complex, k: float) -> dict:
    return _el.rayleigh_cross_sections(alpha_eff, k)


def dipole_t_scalar(alpha: complex, k: float) -> complex:
    return _el.dipole_t_scalar(alpha, k)


def full_wave(*args, **kwargs):
    """Полноволновой конфокальный сфероид (SVM в сфероидальном базисе) — не реализован.

    Требует устойчивых сфероидальных волновых функций с комплексным c (Van Buren
    arXiv:2009.01618 / SMARTIES); scipy поддерживает только вещественный c.
    """
    raise NotImplementedError(
        "Полноволновой сфероид требует библиотеки сфероидальных функций с комплексным c "
        "(scipy pro_*/obl_* — только вещественный c). См. GreenTensor_Theory.tex, раздел «Сфероид»."
    )
