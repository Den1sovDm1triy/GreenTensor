# SPDX-License-Identifier: MIT

"""Независимый аналитический арбитр: рассеяние плоской волны на ОДНОРОДНОМ
бесконечном круговом цилиндре (нормальное падение), Bohren & Huffman гл. 8.

Используется как эталон для проверки многослойного ТФГ-решателя
(green_tensor.cylinder.layered_*). Выведено из первых принципов независимо от
кода библиотеки. Конвенция: e^{-iωt}, исходящая волна ~ H_n^{(1)}, внешняя среда —
вакуум, относительный показатель m = √(εμ).

Поляризации (стандартная конвенция B&H):
  TM = case I  — E параллельно оси цилиндра (E_z),  множитель m/μ;
  TE = case II — H параллельно оси цилиндра (H_z),  множитель m/ε.

Сечения (нормировка на диаметр, на единицу длины), B&H (8.36)-(8.37):
  Q_sca = (2/x) [ |a_0|² + 2 Σ_{n≥1} |a_n|² ],
  Q_ext = (2/x) Re[ a_0 + 2 Σ_{n≥1} a_n ].
"""
from __future__ import annotations

import numpy as np
import scipy.special as sp


def coeff(eps: complex, mu: complex, x: float, n: int, pol: str) -> complex:
    """Коэффициент рассеяния a_n однородного цилиндра (n≥0), pol='TM'|'TE'."""
    m = np.sqrt(eps * mu)
    Jx, Jpx = sp.jv(n, x), sp.jvp(n, x)
    Jm, Jpm = sp.jv(n, m * x), sp.jvp(n, m * x)
    Hx, Hpx = sp.hankel1(n, x), sp.h1vp(n, x)
    if pol == "TM":      # E ∥ ось
        f = m / mu
    elif pol == "TE":    # H ∥ ось
        f = m / eps
    else:
        raise ValueError("pol must be 'TM' or 'TE'")
    num = Jpx * Jm - f * Jpm * Jx
    den = Hpx * Jm - f * Jpm * Hx
    return num / den


def coeffs(eps: complex, mu: complex, x: float, nmax: int, pol: str) -> np.ndarray:
    return np.array([coeff(eps, mu, x, n, pol) for n in range(nmax + 1)], dtype=complex)


def _nmax(x: float) -> int:
    return int(np.ceil(abs(x) + 4.0 * abs(x) ** (1.0 / 3.0) + 2.0)) + 2


def cross_sections(eps: complex, mu: complex, x: float, pol: str,
                   nmax: int | None = None) -> dict:
    """Q_sca, Q_ext, Q_abs однородного цилиндра (нормировка на диаметр)."""
    if nmax is None:
        nmax = _nmax(x)
    a = coeffs(eps, mu, x, nmax, pol)
    q_sca = (2.0 / x) * (abs(a[0]) ** 2 + 2.0 * np.sum(np.abs(a[1:]) ** 2))
    q_ext = (2.0 / x) * np.real(a[0] + 2.0 * np.sum(a[1:]))
    return {"q_sca": float(q_sca), "q_ext": float(q_ext), "q_abs": float(q_ext - q_sca)}


# --------------------------------------------------------------------------- #
# Косое падение (oblique incidence): Wait (1955) / Kavaklıoğlu — независимый
# эталон для однородного диэлектрического цилиндра. Связь E- и H-поляризаций.
# Конвенция e^{-iωt}, H_n^{(1)} исходящая. Нормированные единицы ε0=μ0=1.
#   k_r = sinθ_i, k_z = cosθ_i (в ед. k0); k_l = √(ε_r μ_r − cos²θ_i);
#   ξ = x·sinθ_i (внешний аргумент), η = x·k_l (внутренний), x = k0·a.
# Кросс-связь ∝ n·F·cosθ_i обнуляется при θ_i=90° (нормальное падение).
# --------------------------------------------------------------------------- #
def _oblique_a(eps_r, mu_r, x, theta_i, n, zeta_r):
    """Со-поляризационный коэффициент a_n^ζ (ζ_r = ε_r или μ_r), Kavaklıoğlu ур. (6)."""
    ct = np.cos(theta_i)
    xi = x * np.sin(theta_i)
    kl = np.sqrt(eps_r * mu_r - ct ** 2)
    eta = x * kl
    krl = np.sin(theta_i) / kl                      # k_r / k_l
    Jl, Jpl = sp.jv(n, eta), sp.jvp(n, eta)
    Jr, Jpr = sp.jv(n, xi), sp.jvp(n, xi)
    Hr, Hpr = sp.hankel1(n, xi), sp.h1vp(n, xi)
    num = Jl * Jpr - zeta_r * krl * Jr * Jpl
    den = Jl * Hpr - zeta_r * krl * Hr * Jpl
    return num / den


def _oblique_b(eps_r, mu_r, x, theta_i, n, zeta_r):
    """Кросс-поляризационный коэффициент b_n^ζ, Kavaklıoğlu ур. (7)-(8); ε0=μ0=1."""
    ct = np.cos(theta_i)
    xi = x * np.sin(theta_i)
    kl = np.sqrt(eps_r * mu_r - ct ** 2)
    eta = x * kl
    krl = np.sin(theta_i) / kl
    Jl, Jpl = sp.jv(n, eta), sp.jvp(n, eta)
    Hr, Hpr = sp.hankel1(n, xi), sp.h1vp(n, xi)
    F = (mu_r * eps_r - 1.0) * ct / (mu_r * eps_r - ct ** 2)
    den = Jl * Hpr - zeta_r * krl * Hr * Jpl
    return (Jl * Hr / den) * (1j * n * F / xi)


def oblique_coeffs_TM(eps_r, mu_r, x, theta_i, n):
    """Рассеянные коэффициенты (A_co=E_z, A_cross=H_z) для TM-падения (E∥ось), 1 гармоника.

    Решение связанной 2×2 системы (Kavaklıoğlu ур. 4a,4b для одиночного цилиндра):
        b^μ (A + c·E^i) = −A^H ;   b^ε A^H = A + a^ε·E^i,
    где E^i = sinθ_i (вертикальная поляризация, E0=1, ψ_i=0).
    """
    xi = x * np.sin(theta_i)
    c_n = sp.jv(n, xi) / sp.hankel1(n, xi)
    a_eps = _oblique_a(eps_r, mu_r, x, theta_i, n, eps_r)
    b_mu = _oblique_b(eps_r, mu_r, x, theta_i, n, mu_r)
    b_eps = _oblique_b(eps_r, mu_r, x, theta_i, n, eps_r)
    E_i = np.sin(theta_i)
    A_H = b_mu * (a_eps - c_n) * E_i / (1.0 + b_mu * b_eps)   # кросс-поляризация (H_z)
    A = b_eps * A_H - a_eps * E_i                              # со-поляризация (E_z)
    return A, A_H


if __name__ == "__main__":
    # Самопроверка: для непоглощающего цилиндра Q_abs ≈ 0 (закон сохранения энергии).
    for pol in ("TM", "TE"):
        cs = cross_sections(4.0, 1.0, 2.0, pol)
        assert abs(cs["q_abs"]) < 1e-12, (pol, cs)
        assert cs["q_sca"] > 0
    # Поглощение ⇒ Q_abs > 0.
    cs = cross_sections(4.0 + 1.0j, 1.0, 2.0, "TM")
    assert cs["q_abs"] > 0, cs
    print("analytic_cylinder self-check OK")
