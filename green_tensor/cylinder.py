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
    if mode == "TM":      # E ∥ оси (case I, Bohren & Huffman)
        num = Jmx * Jpx - m * Jx * Jpmx
        den = Jmx * Hpx - m * Hx * Jpmx
    elif mode == "TE":    # H ∥ оси (case II, Bohren & Huffman)
        num = m * Jmx * Jpx - Jx * Jpmx
        den = m * Jmx * Hpx - Hx * Jpmx
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


def _nmax_default(x: float) -> int:
    return int(np.ceil(abs(x) + 4.0 * abs(x) ** (1.0 / 3.0) + 2.0)) + 2


def _layer_basis(n: int, z: complex):
    """(J_n, H_n^{(1)}, J_n', H_n^{(1)'}) при (комплексном) аргументе z.

    Базис {J_n, H_n^{(1)}} устойчив для комплексного z (поглощающие слои); все
    функции scipy (jv/hankel1/jvp/h1vp) поддерживают комплексный аргумент.
    """
    return sp.jv(n, z), sp.hankel1(n, z), sp.jvp(n, z), sp.h1vp(n, z)


def layered_coeff(eps, mu, a_norm, x: float, n: int, mode: str = "TM") -> complex:
    """Коэффициент рассеяния a_n радиально-СЛОИСТОГО цилиндра (нормальное падение).

    Метод эквивалентных линий передачи (ТФГ, Дайлис–Шабунин 2017) при h=0: E- и
    H-линии расцеплены, поэтому каждая поляризация — независимая скалярная
    рекурсия амплитуд (A_i, B_i) поля E_z (TM) или H_z (TE) в базисе {J_n, H_n^{(1)}}.

    eps, mu : проницаемости слоёв изнутри наружу (комплексные допустимы);
    a_norm  : нормированные радиусы границ слоёв (внешний = 1), по возрастанию;
    x       : параметр размера на внешнем радиусе, x = k0·R;
    mode    : 'TM' (E∥ось) или 'TE' (H∥ось).
    """
    eps = [complex(e) for e in eps]
    mu = [complex(m) for m in mu]
    nL = len(eps)
    if len(mu) != nL or len(a_norm) != nL:
        raise ValueError("eps, mu, a_norm must have equal length")
    mrel = [np.sqrt(eps[i] * mu[i]) for i in range(nL)]   # относит. показатель слоя
    # множитель в условии непрерывности тангенциального H (TM) или E (TE):
    fac = [(mrel[i] / mu[i] if mode == "TM" else mrel[i] / eps[i]) for i in range(nL)]

    # ядро регулярно на оси ⇒ только J_n (амплитуда B при H_n^{(1)} равна нулю)
    A, B = 1.0 + 0j, 0.0 + 0j
    for j in range(nL - 1):                # внутренние границы r_j
        u = x * a_norm[j]
        Ji, Hi, Jpi, Hpi = _layer_basis(n, mrel[j] * u)
        Jo, Ho, Jpo, Hpo = _layer_basis(n, mrel[j + 1] * u)
        val = A * Ji + B * Hi                       # непрерывность E_z (TM)/H_z (TE)
        der = fac[j] * (A * Jpi + B * Hpi)          # непрерывность H_φ (TM)/E_φ (TE)
        a21, a22 = fac[j + 1] * Jpo, fac[j + 1] * Hpo
        det = Jo * a22 - Ho * a21
        A = (a22 * val - Ho * der) / det
        B = (-a21 * val + Jo * der) / det

    # внешняя граница u = x (a_norm[-1] = 1): сшивание с вакуумом
    JN, HN, JpN, HpN = _layer_basis(n, mrel[-1] * x)
    G = fac[-1] * (A * JpN + B * HpN) / (A * JN + B * HN)   # поверхностная проводимость
    Jx, Jpx = sp.jv(n, x), sp.jvp(n, x)
    Hx, Hpx = sp.hankel1(n, x), sp.h1vp(n, x)
    return (Jpx - G * Jx) / (Hpx - G * Hx)


def cross_sections_layered(eps, mu, a_norm, x: float, mode: str = "TM",
                           nmax: int | None = None) -> dict:
    """Эффективности слоистого цилиндра (норм. на диаметр): Q_sca, Q_ext, Q_abs.

    При одном слое и μ=1 совпадает с :func:`cross_sections_infinite`.
    """
    if nmax is None:
        nmax = _nmax_default(x)
    a = np.array([layered_coeff(eps, mu, a_norm, x, n, mode) for n in range(nmax + 1)],
                 dtype=complex)
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
