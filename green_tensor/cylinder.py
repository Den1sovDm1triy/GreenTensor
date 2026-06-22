"""cylinder — рассеяние на круговом цилиндре.

БЕСКОНЕЧНЫЙ круговой цилиндр (нормальное падение) решается ТОЧНО разделением
переменных в цилиндрических координатах (2D-ряд Ми, Bohren & Huffman гл. 8):
TM (E вдоль оси) и TE (H вдоль оси) расцеплены, коэффициенты через цилиндрические
Бессель/Ханкель целого порядка. Это каноническая аналитика (одна координатная
поверхность ρ=a).

КОНЕЧНЫЙ цилиндр — несепарабелен (две координатные поверхности + ребро), требует
сшивания мод / EBCM (см. GreenTensor_Theory.tex, раздел «Конечный цилиндр»). Эта
ветка не выдаётся непроверенной заглушкой: фините-решатель поднимает
NotImplementedError, а цилиндроподобные КОНЕЧНЫЕ тела для GMM собираются из
аналитического семейства через decompose.py (стопка сфер) либо аппроксимируются
вытянутым сфероидом (spheroid.py).

Конвенция e^{-iωt}; исходящая волна ~ H_n^{(1)}. x=k·a — параметр размера.
"""
from __future__ import annotations

import numpy as np
import scipy.special as sp


def _coeffs(m: complex, x: float, nmax: int, mode: str):
    """Коэффициенты рассеяния бесконечного цилиндра a_n, n=0..nmax (TM или TE)."""
    n = np.arange(0, nmax + 1)
    Jx, Jpx = sp.jv(n, x), sp.jvp(n, x)
    Jmx, Jpmx = sp.jv(n, m * x), sp.jvp(n, m * x)
    Hx, Hpx = sp.hankel1(n, x), sp.h1vp(n, x)
    if mode == "TM":      # E вдоль оси
        num = m * Jmx * Jpx - Jx * Jpmx
        den = m * Jmx * Hpx - Hx * Jpmx
    elif mode == "TE":    # H вдоль оси
        num = Jmx * Jpx - m * Jx * Jpmx
        den = Jmx * Hpx - m * Hx * Jpmx
    else:
        raise ValueError("mode must be 'TM' or 'TE'")
    return num / den


def cross_sections_infinite(m: complex, x: float, mode: str = "TM", nmax: int | None = None) -> dict:
    """Эффективности бесконечного цилиндра (норм. на диаметр 2a): Q_sca, Q_ext, Q_abs."""
    if nmax is None:
        nmax = int(np.ceil(abs(x) + 4.0 * abs(x) ** (1.0 / 3.0) + 2.0)) + 2
    a = _coeffs(m, x, nmax, mode)
    q_sca = (2.0 / x) * (abs(a[0]) ** 2 + 2.0 * np.sum(np.abs(a[1:]) ** 2))
    q_ext = (2.0 / x) * np.real(a[0] + 2.0 * np.sum(a[1:]))
    return {"q_sca": float(q_sca), "q_ext": float(q_ext), "q_abs": float(q_ext - q_sca)}


def finite(*args, **kwargs):
    """Конечный цилиндр (сшивание мод / EBCM) — не реализован.

    Несепарабелен (ρ=a и z=±L/2 + ребро); требует связанной модальной системы с
    краевым условием Мейкснера. Для GMM конечные цилиндры собираются из сфер
    (decompose.py) или приближаются вытянутым сфероидом (spheroid.py).
    """
    raise NotImplementedError(
        "Конечный цилиндр требует мод-матчинга/EBCM (см. GreenTensor_Theory.tex, "
        "раздел «Конечный цилиндр»). Используйте decompose.py или spheroid.py для GMM."
    )
