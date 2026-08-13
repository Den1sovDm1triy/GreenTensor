# SPDX-License-Identifier: MIT

"""Multilayer radially-anisotropic (radially uniaxial) Mie solver.

Каждый слой l сферы описывается диагональными одноосными по радиусу тензорами

    eps_hat_l = diag(eps_r,l, eps_t,l, eps_t,l),
    mu_hat_l  = diag(mu_r,l,  mu_t,l,  mu_t,l).

Разделение переменных (метод ТФГ, сферическая одноосная анизотропия) оставляет
угловую часть — обычными сферическими гармониками целого порядка n, но радиальные
уравнения приобретают *модифицированный* (в общем случае дробный) угловой порядок
(eq:nu-orders):

    nu_E (nu_E + 1) = (eps_t/eps_r) n (n + 1)      (TM / E-волны),
    nu_H (nu_H + 1) = (mu_t /mu_r ) n (n + 1)      (TE / H-волны),
    nu = -1/2 + sqrt(1/4 + (eps_t/eps_r) n (n+1)).

Трансверсальное волновое число слоя k_t,l = k0 * sqrt(eps_t,l * mu_t,l); радиальные
функции — Риккати-Бессель дробного порядка psi_nu(k_t r), chi_nu(k_t r), xi_nu(k_t r).

Сшивание слоёв — тот же рекуррентный пересчёт направленных импедансов Z (TM) и
адмитансов Y (TE), что и в изотропном многослойном ядре 01_sphere.py
(RCSCalculator.calculate_impedances), но с порядком nu_E / nu_H вместо целого n и
трансверсальными eps_t, mu_t в нормировках (eq:Z-recur-aniso):

    Z_0 = (eta_0/eta_1) * psi'_{nu_E,0}(m_0 v_0) / psi_{nu_E,0}(m_0 v_0),
    Z_h = (eta_h/eta_{h+1}) * (Cpr_h + Z_{h-1} Spr_h) / (C_h + Z_{h-1} S_h),
    eta_l = sqrt(mu_t,l / eps_t,l)   (TM),   Y-ветвь дуальна: sqrt(eps_t,l/mu_t,l).

C_h, S_h, Cpr_h, Spr_h — те же вронскианоподобные передаточные комбинации оболочки h
(порядок nu_h, волновое число m_h = sqrt(eps_t,h mu_t,h)), что в 01_sphere.py, но с
дробным порядком. Внешняя среда — вакуум (eta = 1), функции целого порядка n.

Замыкание на внешней границе (x = k0 R, вакуум, порядок n):

    a_n = (Z_out psi_n(x) - psi'_n(x)) / (Z_out xi_n(x) - xi'_n(x)),
    b_n = (Y_out psi_n(x) - psi'_n(x)) / (Y_out xi_n(x) - xi'_n(x)).

В пределе eps_r=eps_t, mu_r=mu_t (nu -> n) рекурсия тождественно совпадает с
магнитным многослойным ядром 01_sphere.py; при L=1 и mu_t=1 — с одиночной
анизотропной сферой (_dissSource/.../anisotropic_mie.py, формула Geng et al. 2004).

Конвенция: e^{-i omega t}, исходящая волна ~ h_n^{(1)}, xi_n(z) = z h_n^{(1)}(z).
Коэффициенты a_n, b_n возвращаются в этой конвенции (как в anisotropic_mie.py);
каноническое ядро 01_sphere.py даёт их комплексно-сопряжёнными (Mn=a_n*, Nn=b_n*),
что не влияет на Q_sca, Q_ext, Q_back.

Ссылка на аналогичные формулы одиночной одноосной сферы:
  Geng Y.-L., Wu X.-B., Li L.-W., Guan B.-R. Mie scattering by a uniaxial
  anisotropic sphere // Phys. Rev. E. 2004. Vol. 70. 056609.
"""
from __future__ import annotations

import numpy as np
import scipy.special as sp

__all__ = [
    "angular_orders",
    "mie_multilayer_anisotropic",
    "cross_sections_anisotropic",
    "scattering_amplitudes_anisotropic",
    "bistatic_pattern_anisotropic",
]


# --------------------------------------------------------------------------- #
# Riccati functions of arbitrary (fractional) order for complex argument
# --------------------------------------------------------------------------- #
def _riccati_nu(nu: float, z):
    """Риккати-Бессель дробного порядка nu для (комплексного) аргумента z.

    Возвращает (psi, psi', chi, chi', xi, xi'), где
        psi_nu(z) = z j_nu(z) = sqrt(pi z / 2) J_{nu+1/2}(z),
        chi_nu(z) = z y_nu(z) = sqrt(pi z / 2) Y_{nu+1/2}(z),
        xi_nu(z)  = z h_nu^{(1)}(z) = psi_nu + i chi_nu,
        f'_nu(z)  = sqrt(pi z / 2) F_{nu-1/2}(z) - (nu/z) f_nu(z).

    Формула производной — обобщение _riccati из tests/analytic_mie.py на дробный
    порядок (DLMF 10.51: F'_{nu}(z) = F_{nu-1}(z) - ((nu+1)/z) F_nu(z), для
    сферических функций даёт psi'_nu = sqrt(pi z/2) J_{nu-1/2} - (nu/z) psi_nu).
    Через jv/yv/hankel1 полуцелого порядка — корректно для комплексного z
    (поглощающие слои) и любого вещественного nu.
    """
    z = np.asarray(z, dtype=complex)
    s = np.sqrt(np.pi * z / 2.0)
    a = nu + 0.5
    am = nu - 0.5
    psi = s * sp.jv(a, z)
    psi_p = s * sp.jv(am, z) - (nu / z) * psi
    chi = s * sp.yv(a, z)
    chi_p = s * sp.yv(am, z) - (nu / z) * chi
    xi = s * sp.hankel1(a, z)
    xi_p = s * sp.hankel1(am, z) - (nu / z) * xi
    return psi, psi_p, chi, chi_p, xi, xi_p


# --------------------------------------------------------------------------- #
# Modified angular orders
# --------------------------------------------------------------------------- #
def angular_orders(n: int, eps_r, eps_t, mu_r, mu_t) -> tuple[float, float]:
    """Модифицированные угловые порядки (nu_E, nu_H) слоя для мультиполя n (eq:nu-orders).

        nu_E = -1/2 + sqrt(1/4 + (eps_t/eps_r) n(n+1))   — TM/E,
        nu_H = -1/2 + sqrt(1/4 + (mu_t /mu_r ) n(n+1))   — TE/H.

    Порядки вещественны (собственные значения сферического лапласиана): зависят
    только от вещественных положительных отношений eps_t/eps_r и mu_t/mu_r.
    При eps_r=eps_t, mu_r=mu_t дают nu_E=nu_H=n (изотропный предел).
    """
    rE = float(np.real(eps_t / eps_r))
    rH = float(np.real(mu_t / mu_r))
    nn = n * (n + 1.0)
    nu_E = -0.5 + np.sqrt(0.25 + rE * nn)
    nu_H = -0.5 + np.sqrt(0.25 + rH * nn)
    return float(nu_E), float(nu_H)


# --------------------------------------------------------------------------- #
# Helpers: broadcast material inputs to per-layer arrays
# --------------------------------------------------------------------------- #
def _as_layer_array(value, L: int, name: str) -> np.ndarray:
    arr = np.atleast_1d(np.asarray(value, dtype=complex))
    if arr.size == 1:
        arr = np.full(L, arr[0], dtype=complex)
    if arr.size != L:
        raise ValueError(f"{name}: ожидалось {L} слоёв, получено {arr.size}")
    return arr


# --------------------------------------------------------------------------- #
# Core solver — multilayer anisotropic Mie coefficients
# --------------------------------------------------------------------------- #
def mie_multilayer_anisotropic(
    k0: float,
    a_norm,
    R: float,
    eps_r,
    eps_t,
    mu_r,
    mu_t,
    nmax: int | None = None,
):
    """Коэффициенты рассеяния (a_n, b_n) многослойной радиально-одноосной сферы.

    Parameters
    ----------
    k0 : float
        Волновое число внешней среды (вакуум). Параметр размера x = k0 * R.
    a_norm : sequence
        Нормированные радиусы границ слоёв ИЗНУТРИ НАРУЖУ, по возрастанию,
        внешняя граница = 1.0. Длина = число слоёв L.
    R : float
        Внешний радиус сферы.
    eps_r, eps_t, mu_r, mu_t : scalar | sequence(len=L)
        Диагональные компоненты тензоров по слоям (относительно вакуума).
        Скаляр разворачивается на все слои.
    nmax : int, optional
        Число членов ряда; по умолчанию правило Уискомба x + 4 x^{1/3} + 2.

    Returns
    -------
    a, b : np.ndarray (complex, длина nmax)
        Коэффициенты Ми (TM=a_n, TE=b_n), индексация n=1..nmax (a[0]=a_1).
        Конвенция e^{-i omega t}, исходящая волна h_n^{(1)}.
    """
    a_norm = np.asarray(a_norm, dtype=float).ravel()
    L = a_norm.size
    if L < 1:
        raise ValueError("нужен хотя бы один слой")
    if np.any(np.diff(a_norm) <= 0):
        raise ValueError("a_norm должен строго возрастать (изнутри наружу)")
    if not np.isclose(a_norm[-1], 1.0):
        raise ValueError("внешняя граница a_norm[-1] должна быть 1.0")

    eps_r = _as_layer_array(eps_r, L, "eps_r")
    eps_t = _as_layer_array(eps_t, L, "eps_t")
    mu_r = _as_layer_array(mu_r, L, "mu_r")
    mu_t = _as_layer_array(mu_t, L, "mu_t")

    x = float(k0) * float(R)                       # внешний параметр размера
    if nmax is None:
        nmax = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
    nmax = max(int(nmax), 1)

    # Трансверсальный показатель преломления и волновые импедансы/адмитансы слоёв.
    m = np.sqrt(eps_t * mu_t)                       # k_t,l = k0 * m_l
    eta = np.sqrt(mu_t / eps_t)                     # TM импеданс слоя (gamma_e)
    Yad = np.sqrt(eps_t / mu_t)                     # TE адмитанс слоя (gamma_h)
    # Внешняя среда (вакуум): eta = Y = 1.
    eta_ext = np.append(eta, 1.0)
    Yad_ext = np.append(Yad, 1.0)

    # Размерные аргументы границ: во внешней (вакуумной) шкале s_l = k0 R a_norm[l];
    # внутри слоя l аргумент на его границах = m_l * s.
    s = x * a_norm                                  # s_l = k0 R a_norm[l]

    a = np.zeros(nmax, dtype=complex)
    b = np.zeros(nmax, dtype=complex)

    for n in range(1, nmax + 1):
        # Порядки nu_E, nu_H по слоям для данного мультиполя n.
        nuE = np.empty(L)
        nuH = np.empty(L)
        for l in range(L):
            nuE[l], nuH[l] = angular_orders(n, eps_r[l], eps_t[l], mu_r[l], mu_t[l])

        Z = _impedance_recursion(nuE, m, s, eta_ext)   # TM
        Y = _impedance_recursion(nuH, m, s, Yad_ext)   # TE

        # Замыкание на вакуум (порядок n, аргумент x).
        psi_x, psi_xp, _, _, xi_x, xi_xp = _riccati_nu(float(n), x)
        psi_x = psi_x.item(); psi_xp = psi_xp.item()
        xi_x = xi_x.item(); xi_xp = xi_xp.item()

        a[n - 1] = (Z * psi_x - psi_xp) / (Z * xi_x - xi_xp)
        b[n - 1] = (Y * psi_x - psi_xp) / (Y * xi_x - xi_xp)

    return a, b


def _impedance_recursion(nu_layers: np.ndarray, m: np.ndarray, s: np.ndarray,
                         gamma_ext: np.ndarray) -> complex:
    """Направленная рекурсия импеданса/адмитанса наружу через L слоёв (eq:Z-recur-aniso).

    nu_layers : порядки nu слоёв (nu_E для TM, nu_H для TE);
    m         : трансверсальные показатели преломления слоёв;
    s         : размерные аргументы границ во внешней шкале (s_l = k0 R a_norm[l]);
    gamma_ext : импедансы/адмитансы слоёв + внешняя среда (длина L+1, gamma_ext[L]=1).

    Возвращает поверхностный импеданс/адмитанс на внешней границе.
    """
    L = nu_layers.size

    # Оболочка 0 (ядро, регулярно в центре): только psi.
    beta0 = m[0] * s[0]
    psi0, psi0_p, _, _, _, _ = _riccati_nu(nu_layers[0], beta0)
    Z = (gamma_ext[0] / gamma_ext[1]) * (psi0_p / psi0)
    Z = complex(Z)

    # Оболочки 1..L-1: передаточные комбинации Вронскиана + скачок импеданса.
    for h in range(1, L):
        alpha = m[h] * s[h - 1]     # внутренняя граница оболочки h
        beta = m[h] * s[h]          # внешняя граница оболочки h
        nu = nu_layers[h]
        pa, pa_p, ca, ca_p, _, _ = _riccati_nu(nu, alpha)
        pb, pb_p, cb, cb_p, _, _ = _riccati_nu(nu, beta)
        pa = pa.item(); pa_p = pa_p.item(); ca = ca.item(); ca_p = ca_p.item()
        pb = pb.item(); pb_p = pb_p.item(); cb = cb.item(); cb_p = cb_p.item()

        C = pb * ca_p - cb * pa_p
        Cpr = pb_p * ca_p - cb_p * pa_p
        S = cb * pa - pb * ca
        Spr = cb_p * pa - pb_p * ca

        Z = (gamma_ext[h] / gamma_ext[h + 1]) * (Cpr + Z * Spr) / (C + Z * S)

    return complex(Z)


# --------------------------------------------------------------------------- #
# Cross sections
# --------------------------------------------------------------------------- #
def cross_sections_anisotropic(
    k0: float,
    a_norm,
    R: float,
    eps_r,
    eps_t,
    mu_r,
    mu_t,
    nmax: int | None = None,
) -> dict:
    """Эффективности Q_sca, Q_ext, Q_abs, Q_back (нормировка на pi R^2).

        Q_sca  = (2/x^2) sum (2n+1) (|a_n|^2 + |b_n|^2),
        Q_ext  = (2/x^2) sum (2n+1) Re(a_n + b_n),
        Q_back = (1/x^2) |sum (2n+1) (-1)^n (a_n - b_n)|^2,
        x = k0 R.
    """
    a, b = mie_multilayer_anisotropic(k0, a_norm, R, eps_r, eps_t, mu_r, mu_t, nmax)
    x = float(k0) * float(R)
    n = np.arange(1, len(a) + 1)
    q_sca = float((2.0 / x**2) * np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2)))
    q_ext = float((2.0 / x**2) * np.sum((2 * n + 1) * np.real(a + b)))
    q_back = float((1.0 / x**2) * np.abs(np.sum((2 * n + 1) * ((-1) ** n) * (a - b))) ** 2)
    return {"Q_sca": q_sca, "Q_ext": q_ext, "Q_abs": q_ext - q_sca, "Q_back": q_back}


# --------------------------------------------------------------------------- #
# Angular pattern (bistatic scattering diagram)
# --------------------------------------------------------------------------- #
def _pi_tau(nmax: int, theta: np.ndarray):
    """Угловые функции pi_n, tau_n (восходящая рекуррентность Bohren & Huffman)."""
    mu = np.cos(np.asarray(theta, dtype=float))
    T = mu.shape[0]
    pin = np.zeros((nmax, T))
    taun = np.zeros((nmax, T))
    pin[0] = 1.0
    taun[0] = mu
    if nmax >= 2:
        pin[1] = 3.0 * mu
        taun[1] = 2.0 * mu * pin[1] - 3.0 * pin[0]
    for n in range(3, nmax + 1):
        i = n - 1
        pin[i] = ((2 * n - 1) / (n - 1)) * mu * pin[i - 1] - (n / (n - 1)) * pin[i - 2]
        taun[i] = n * mu * pin[i] - (n + 1) * pin[i - 1]
    return pin, taun


def scattering_amplitudes_anisotropic(
    k0: float, a_norm, R: float, eps_r, eps_t, mu_r, mu_t,
    theta: np.ndarray, nmax: int | None = None,
):
    """Комплексные амплитуды рассеяния S1(theta), S2(theta) (линейная поляризация).

        S1 = sum (2n+1)/(n(n+1)) (a_n pi_n + b_n tau_n),
        S2 = sum (2n+1)/(n(n+1)) (a_n tau_n + b_n pi_n).
    """
    a, b = mie_multilayer_anisotropic(k0, a_norm, R, eps_r, eps_t, mu_r, mu_t, nmax)
    nmax = len(a)
    pin, taun = _pi_tau(nmax, theta)
    n = np.arange(1, nmax + 1)
    w = ((2 * n + 1) / (n * (n + 1)))[:, None]
    ac, bc = a[:, None], b[:, None]
    S1 = np.sum(w * (ac * pin + bc * taun), axis=0)
    S2 = np.sum(w * (ac * taun + bc * pin), axis=0)
    return S1, S2


def bistatic_pattern_anisotropic(
    k0: float, a_norm, R: float, eps_r, eps_t, mu_r, mu_t,
    theta: np.ndarray, nmax: int | None = None,
) -> np.ndarray:
    """Бистатическая диаграмма рассеяния (|S1|^2 + |S2|^2)/x^2, неполяризованная."""
    S1, S2 = scattering_amplitudes_anisotropic(
        k0, a_norm, R, eps_r, eps_t, mu_r, mu_t, theta, nmax)
    x = float(k0) * float(R)
    return (np.abs(S1) ** 2 + np.abs(S2) ** 2) / x**2
