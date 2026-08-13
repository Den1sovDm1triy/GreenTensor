# SPDX-License-Identifier: MIT

"""Независимый аналитический эталон Ми для однородной сферы (Bohren & Huffman).

Используется как «арбитр корректности» для расчётов GreenTensor: коэффициенты
a_n, b_n, сечения Q_sca / Q_ext / Q_abs и PEC-предел выводятся из замкнутых формул
через цилиндрические функции scipy и не зависят от кода библиотеки.

Конвенция: e^{-iωt}, ξ_n(x) = x·h_n^{(1)}(x) (исходящая волна).
Поглощающая среда: Im(m) > 0, т.е. Im(eps) > 0  (m = sqrt(eps), eps = eps' + i·eps'').

Спецфункции вычисляются через jv/yv/hankel1 (полуцелый порядок) — они принимают
КОМПЛЕКСНЫЙ аргумент, поэтому эталон корректен для поглощающих сред и металлов
(в отличие от scipy.special.spherical_jn, работающего только с вещественным z).
"""
from __future__ import annotations

import numpy as np
import scipy.special as sp

# np.trapz переименован в np.trapezoid в numpy 2.0; поддерживаем обе версии
trapezoid = getattr(np, "trapezoid", getattr(np, "trapz"))


def wiscombe_nmax(x: float) -> int:
    """Число членов ряда Ми по правилу Уискомба."""
    x = float(np.real(x)) if np.iscomplexobj(x) else float(x)
    return int(np.ceil(abs(x) + 4.0 * abs(x) ** (1.0 / 3.0) + 2.0))


def _riccati(n: np.ndarray, z: complex):
    """Риккати-Бессель ψ_n, ψ'_n и ξ_n, ξ'_n для (комплексного) аргумента z.

    ψ_n(z) = z·j_n(z) = sqrt(πz/2)·J_{n+1/2}(z)
    ξ_n(z) = z·h_n^{(1)}(z) = sqrt(πz/2)·H^{(1)}_{n+1/2}(z)
    ψ'_n(z) = sqrt(πz/2)·J_{n-1/2}(z) − (n/z)·ψ_n(z)   (то же для ξ через H^{(1)})
    """
    s = np.sqrt(np.pi * z / 2.0)
    nu = n + 0.5
    num = n - 0.5
    psi = s * sp.jv(nu, z)
    psi_p = s * sp.jv(num, z) - (n / z) * psi
    xi = s * sp.hankel1(nu, z)
    xi_p = s * sp.hankel1(num, z) - (n / z) * xi
    return psi, psi_p, xi, xi_p


def mie_ab(m: complex, x: float, nmax: int | None = None):
    """Коэффициенты Ми a_n, b_n однородной сферы (m может быть комплексным)."""
    if nmax is None:
        nmax = wiscombe_nmax(x)
    n = np.arange(1, nmax + 1)
    psi_x, psi_xp, xi_x, xi_xp = _riccati(n, x)
    psi_mx, psi_mxp, _, _ = _riccati(n, m * x)
    a = (m * psi_mx * psi_xp - psi_x * psi_mxp) / (m * psi_mx * xi_xp - xi_x * psi_mxp)
    b = (psi_mx * psi_xp - m * psi_x * psi_mxp) / (psi_mx * xi_xp - m * xi_x * psi_mxp)
    return n, a, b


def mie_ab_eps_mu(eps: complex, mu: complex, x: float, nmax: int | None = None):
    """Коэффициенты Ми a_n, b_n для однородной магнитодиэлектрической сферы.

    Это независимая закрытая формула для относительных eps и mu в вакууме. При
    mu=1 она сводится к ``mie_ab(sqrt(eps), x, nmax)``. Запись следует
    формулам Керкера для магнитных сфер:
    m=sqrt(eps*mu), m_tilde=sqrt(eps/mu).
    """
    if nmax is None:
        nmax = wiscombe_nmax(x)
    n = np.arange(1, nmax + 1)
    eps = complex(eps)
    mu = complex(mu)
    m = np.sqrt(eps * mu)
    m_tilde = np.sqrt(eps / mu)
    psi_x, psi_xp, xi_x, xi_xp = _riccati(n, x)
    psi_mx, psi_mxp, _, _ = _riccati(n, m * x)
    a = (m_tilde * psi_mx * psi_xp - psi_x * psi_mxp) / (
        m_tilde * psi_mx * xi_xp - xi_x * psi_mxp
    )
    b = (psi_mx * psi_xp - m_tilde * psi_x * psi_mxp) / (
        psi_mx * xi_xp - m_tilde * xi_x * psi_mxp
    )
    return n, a, b


def mie_pec(x: float, nmax: int | None = None):
    """Коэффициенты Ми для идеального проводника (PEC), предел m → ∞."""
    if nmax is None:
        nmax = wiscombe_nmax(x)
    n = np.arange(1, nmax + 1)
    psi, psi_p, xi, xi_p = _riccati(n, x)
    a = psi_p / xi_p
    b = psi / xi
    return n, a, b


def _q_sca(n, a, b, x):
    return float((2.0 / x**2) * np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2)))


def _q_ext(n, a, b, x):
    return float((2.0 / x**2) * np.sum((2 * n + 1) * np.real(a + b)))


def q_sca(m: complex, x: float, nmax: int | None = None) -> float:
    return _q_sca(*mie_ab(m, x, nmax), x)


def q_ext(m: complex, x: float, nmax: int | None = None) -> float:
    return _q_ext(*mie_ab(m, x, nmax), x)


def q_abs(m: complex, x: float, nmax: int | None = None) -> float:
    n, a, b = mie_ab(m, x, nmax)
    return _q_ext(n, a, b, x) - _q_sca(n, a, b, x)


def q_sca_pec(x: float, nmax: int | None = None) -> float:
    return _q_sca(*mie_pec(x, nmax), x)


def q_sca_rayleigh(m: complex, x: float) -> float:
    """Релеевский предел Q_sca при x << 1 (для самопроверки эталона)."""
    return float((8.0 / 3.0) * x**4 * np.abs((m**2 - 1) / (m**2 + 2)) ** 2)


def q_abs_rayleigh(m: complex, x: float) -> float:
    """Релеевский предел Q_abs при x << 1 (поглощение)."""
    return float(4.0 * x * np.imag((m**2 - 1) / (m**2 + 2)))


# --------------------------------------------------------------------------- #
# Угловые функции и амплитуды рассеяния (НЕЗАВИСИМО: рекуррентность B&H, не lpmv)
# --------------------------------------------------------------------------- #
def mie_pi_tau(nmax: int, theta: np.ndarray):
    """Угловые функции π_n, τ_n по восходящей рекуррентности (Bohren & Huffman).

    π_1=1; π_n = ((2n-1)/(n-1))·μ·π_{n-1} − (n/(n-1))·π_{n-2};
    τ_n = n·μ·π_n − (n+1)·π_{n-1};   μ = cosθ.
    Возвращает массивы [nmax, len(theta)].
    """
    mu = np.cos(np.asarray(theta, dtype=float))
    T = mu.shape[0]
    pin = np.zeros((nmax, T))
    taun = np.zeros((nmax, T))
    pin[0] = 1.0
    taun[0] = mu  # 1·μ·1 − 2·π_0
    if nmax >= 2:
        pin[1] = 3.0 * mu  # n=2
        taun[1] = 2.0 * mu * pin[1] - 3.0 * pin[0]
    for n in range(3, nmax + 1):
        i = n - 1
        pin[i] = ((2 * n - 1) / (n - 1)) * mu * pin[i - 1] - (n / (n - 1)) * pin[i - 2]
        taun[i] = n * mu * pin[i] - (n + 1) * pin[i - 1]
    return pin, taun


def mie_amplitudes(m: complex, x: float, theta: np.ndarray, nmax: int | None = None):
    """Линейные амплитуды рассеяния S1 (перпендикулярная), S2 (параллельная)."""
    if nmax is None:
        nmax = wiscombe_nmax(x)
    n, a, b = mie_ab(m, x, nmax)
    pin, taun = mie_pi_tau(nmax, theta)
    wl = ((2 * n + 1) / (n * (n + 1)))[:, None]
    a_c, b_c = a[:, None], b[:, None]
    S1 = np.sum(wl * (a_c * pin + b_c * taun), axis=0)
    S2 = np.sum(wl * (a_c * taun + b_c * pin), axis=0)
    return S1, S2


def mie_helicity(m: complex, x: float, theta: np.ndarray, nmax: int | None = None):
    """Круговые амплитуды: S_co=(S1+S2)/2 (со-пол.), S_cross=(S1−S2)/2 (кросс-пол.)."""
    S1, S2 = mie_amplitudes(m, x, theta, nmax)
    return 0.5 * (S1 + S2), 0.5 * (S1 - S2)


if __name__ == "__main__":
    # 1) Релеевский предел рассеяния (вещественный m)
    m = 1.5
    for x in (0.01, 0.02, 0.05):
        rel = abs(q_sca(m, x, 4) - q_sca_rayleigh(m, x)) / q_sca_rayleigh(m, x)
        assert rel < 0.03, f"Q_sca не сходится к Рэлею: {rel}"
    # 2) Поглощение: Q_abs > 0 при Im(m) > 0 и согласуется с релеевским пределом
    m = 1.5 + 0.01j
    for x in (0.01, 0.02, 0.05):
        qa = q_abs(m, x, 4)
        ref = q_abs_rayleigh(m, x)
        assert qa > 0, "Q_abs должно быть > 0 для поглощающей среды (Im m > 0)"
        assert abs(qa - ref) / abs(ref) < 0.05, f"Q_abs не сходится к Рэлею: {qa} vs {ref}"
    # 3) PEC: при больших x Q_ext → 2 (парадокс экстинкции); Q_abs(PEC) = 0
    n, a, b = mie_pec(20.0)
    qext = _q_ext(n, a, b, 20.0)
    qsca = _q_sca(n, a, b, 20.0)
    assert abs(qext - 2.0) < 0.05, f"PEC Q_ext при x=20 должно быть ~2, получено {qext}"
    assert abs(qext - qsca) < 1e-6, "PEC лосслесс: Q_ext должно равняться Q_sca"
    # 4) Амплитуды/угловые функции: forward S1≈S2 (кросс-пол=0); ∫(|S1|²+|S2|²)sinθ = x²·Q_sca
    th = np.linspace(1e-4, np.pi, 4000)
    S1, S2 = mie_amplitudes(1.5, 3.0, th)
    integ = trapezoid((np.abs(S1) ** 2 + np.abs(S2) ** 2) * np.sin(th), th)
    assert abs(integ - 9.0 * q_sca(1.5, 3.0)) < 1e-3 * 9.0 * q_sca(1.5, 3.0), "∫|S|² != x²·Qsca"
    Sf1, Sf2 = mie_amplitudes(1.5, 3.0, np.array([1e-3]))
    assert abs(Sf1[0] - Sf2[0]) / abs(Sf1[0]) < 1e-3, "forward: S1 должно равняться S2"
    print("OK: эталон согласован — Рэлей, PEC, амплитуды S1/S2 (forward и интеграл Q_sca).")
