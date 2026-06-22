"""vswf — сферические волновые функции и теоремы сложения (фундамент GMM).

Общий вычислительный слой для сборки сложной геометрии (см. GreenTensor_Theory.tex,
раздел «Суперпозиция T-матриц»). На этом этапе реализован СКАЛЯРНЫЙ слой:

  * символы Вигнера 3j (формула Рака, точные факториалы);
  * нормированные сферические гармоники Y_n^m (с фазой Кондона–Шортли, через lpmv);
  * коэффициенты Гаунта  G(l1,m1;l2,m2;l3,m3) = ∫ Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} dΩ;
  * скалярные регулярные/исходящие волны  ψ^{(1,3)}_{nm} = {j_n|h_n^{(1)}}(kr) Y_n^m;
  * скалярные коэффициенты трансляции (теорема сложения) для регулярных волн:

        Rg ψ_{νμ}(r − d) = Σ_{n,m} α^{νμ}_{nm}(d) · Rg ψ_{nm}(r).

Вывод α — из разложения Рэлея плоской волны; предел d→0 даёт α = δ (проверяется тестом).
Конвенция e^{-iωt}; исходящая волна ~ h_n^{(1)}. Векторные M,N и их трансляции (A,B) —
следующий слой поверх этого (проверенного) фундамента.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
import scipy.special as sp


# --------------------------------------------------------------------------- #
# Символ Вигнера 3j (формула Рака)
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=None)
def wigner_3j(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    """Символ Вигнера 3j для целых аргументов (формула Рака, точные факториалы)."""
    if m1 + m2 + m3 != 0:
        return 0.0
    if not (abs(j1 - j2) <= j3 <= j1 + j2):
        return 0.0
    if any(abs(m) > j for m, j in ((m1, j1), (m2, j2), (m3, j3))):
        return 0.0
    f = math.factorial
    tri = f(j1 + j2 - j3) * f(j1 - j2 + j3) * f(-j1 + j2 + j3) / f(j1 + j2 + j3 + 1)
    pref = tri * (f(j1 + m1) * f(j1 - m1) * f(j2 + m2) * f(j2 - m2)
                  * f(j3 + m3) * f(j3 - m3))
    kmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    kmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    s = 0.0
    for k in range(kmin, kmax + 1):
        denom = (f(k) * f(j1 + j2 - j3 - k) * f(j1 - m1 - k) * f(j2 + m2 - k)
                 * f(j3 - j2 + m1 + k) * f(j3 - j1 - m2 + k))
        s += (-1) ** k / denom
    return (-1) ** (j1 - j2 - m3) * math.sqrt(pref) * s


def gaunt(l1: int, m1: int, l2: int, m2: int, l3: int, m3: int) -> float:
    """Коэффициент Гаунта G = ∫ Y_{l1}^{m1} Y_{l2}^{m2} Y_{l3}^{m3} dΩ (вещественный)."""
    if m1 + m2 + m3 != 0:
        return 0.0
    if (l1 + l2 + l3) % 2 != 0:
        return 0.0
    if not (abs(l1 - l2) <= l3 <= l1 + l2):
        return 0.0
    pre = math.sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4.0 * math.pi))
    return pre * wigner_3j(l1, l2, l3, 0, 0, 0) * wigner_3j(l1, l2, l3, m1, m2, m3)


# --------------------------------------------------------------------------- #
# Сферические гармоники (фаза Кондона–Шортли, как в scipy.lpmv)
# --------------------------------------------------------------------------- #
def ynm(n: int, m: int, theta, phi):
    """Нормированная сферическая гармоника Y_n^m(θ,φ); поддерживает массивы θ,φ."""
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    am = abs(m)
    norm = math.sqrt((2 * n + 1) / (4.0 * math.pi)
                     * math.factorial(n - am) / math.factorial(n + am))
    P = sp.lpmv(am, n, np.cos(theta))  # включает фазу Кондона–Шортли
    Y = norm * P * np.exp(1j * am * phi)
    if m < 0:
        Y = (-1) ** am * np.conj(Y)
    return Y


# --------------------------------------------------------------------------- #
# Скалярные волновые функции
# --------------------------------------------------------------------------- #
def _to_spherical(xyz):
    """(...,3) -> (r, theta, phi). theta∈[0,π], phi∈[0,2π)."""
    xyz = np.asarray(xyz, dtype=float)
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    r = np.sqrt(x * x + y * y + z * z)
    with np.errstate(invalid="ignore", divide="ignore"):
        theta = np.arccos(np.where(r > 0, z / r, 1.0))
    phi = np.arctan2(y, x)
    return r, theta, phi


def rg_scalar(n: int, m: int, k: float, xyz):
    """Регулярная скалярная волна Rg ψ_{nm} = j_n(kr) Y_n^m(θ,φ)."""
    r, theta, phi = _to_spherical(xyz)
    return sp.spherical_jn(n, k * r) * ynm(n, m, theta, phi)


def out_scalar(n: int, m: int, k: float, xyz):
    """Исходящая скалярная волна ψ^{(3)}_{nm} = h_n^{(1)}(kr) Y_n^m(θ,φ)."""
    r, theta, phi = _to_spherical(xyz)
    h = sp.spherical_jn(n, k * r) + 1j * sp.spherical_yn(n, k * r)
    return h * ynm(n, m, theta, phi)


# --------------------------------------------------------------------------- #
# Скалярная теорема сложения (трансляция регулярных волн)
# --------------------------------------------------------------------------- #
def alpha_scalar(nu: int, mu: int, n: int, m: int, d, k: float) -> complex:
    """Коэффициент трансляции α^{νμ}_{nm}(d): Rg ψ_{νμ}(r−d)=Σ α Rg ψ_{nm}(r)."""
    dr, dth, dph = _to_spherical(d)
    dr = float(dr)
    acc = 0.0 + 0.0j
    s = m - mu
    for p in range(abs(n - nu), n + nu + 1):
        if (n + nu + p) % 2 != 0:      # G(.,.,.;0,0,0) = 0 при нечётной сумме
            continue
        if abs(s) > p:
            continue
        G = gaunt(n, -m, nu, mu, p, s)
        if G == 0.0:
            continue
        Yp = ynm(p, s, dth, dph)
        acc += (-1j) ** p * sp.spherical_jn(p, k * dr) * np.conj(Yp) * G
    return 4.0 * math.pi * (1j) ** (n - nu) * (-1) ** m * acc


# --------------------------------------------------------------------------- #
# Векторные сферические волновые функции M_{nm}, N_{nm}
# --------------------------------------------------------------------------- #
def _dtheta_ynm(n: int, m: int, theta: float, phi: float) -> complex:
    """∂_θ Y_n^m в точке (θ,φ). Через DLMF 14.10.5; θ не на полюсе."""
    if m < 0:
        am = -m
        return (-1) ** am * np.conj(_dtheta_ynm(n, am, theta, phi))
    ct, st = math.cos(theta), math.sin(theta)
    Pn = sp.lpmv(m, n, ct)
    Pn1 = sp.lpmv(m, n - 1, ct) if (n - 1) >= m else 0.0
    dP = -((n + m) * Pn1 - n * ct * Pn) / st           # dP_n^m/dθ
    norm = math.sqrt((2 * n + 1) / (4.0 * math.pi)
                     * math.factorial(n - m) / math.factorial(n + m))
    return norm * dP * np.exp(1j * m * phi)


def _sph_to_cart_vec(Ar, At, Ap, theta: float, phi: float):
    """Вектор (A_r, A_θ, A_φ) в точке (θ,φ) -> декартовы (Ax, Ay, Az)."""
    ct, st, cp, sp_ = math.cos(theta), math.sin(theta), math.cos(phi), math.sin(phi)
    rhat = np.array([st * cp, st * sp_, ct])
    that = np.array([ct * cp, ct * sp_, -st])
    phat = np.array([-sp_, cp, 0.0])
    return Ar * rhat + At * that + Ap * phat


def _zfuncs(n: int, rho: float, kind: str):
    """z_n(ρ), u_n'(ρ)=d(ρ z_n)/dρ для регулярной (j_n) или исходящей (h_n^{(1)}) волны."""
    jn = sp.spherical_jn(n, rho)
    jnp = sp.spherical_jn(n, rho, derivative=True)
    if kind == "reg":
        z, zp = jn, jnp
    elif kind == "out":
        yn = sp.spherical_yn(n, rho)
        ynp = sp.spherical_yn(n, rho, derivative=True)
        z, zp = jn + 1j * yn, jnp + 1j * ynp
    else:
        raise ValueError("kind must be 'reg' or 'out'")
    return z, z + rho * zp     # u_n' = z_n + ρ z_n'


def vsw_M(n: int, m: int, k: float, xyz, kind: str = "reg"):
    """Векторная сферическая волновая функция M_{nm} (декартов вектор) в точке xyz."""
    r, theta, phi = (float(v) for v in _to_spherical(xyz))
    rho = k * r
    z, _ = _zfuncs(n, rho, kind)
    Y = complex(ynm(n, m, theta, phi))
    dY = _dtheta_ynm(n, m, theta, phi)
    st = math.sin(theta)
    M_t = z * (1j * m / st) * Y
    M_p = -z * dY
    return _sph_to_cart_vec(0.0, M_t, M_p, theta, phi)


def vsw_N(n: int, m: int, k: float, xyz, kind: str = "reg"):
    """Векторная сферическая волновая функция N_{nm} (декартов вектор) в точке xyz."""
    r, theta, phi = (float(v) for v in _to_spherical(xyz))
    rho = k * r
    z, up = _zfuncs(n, rho, kind)
    Y = complex(ynm(n, m, theta, phi))
    dY = _dtheta_ynm(n, m, theta, phi)
    st = math.sin(theta)
    N_r = n * (n + 1) * (z / rho) * Y
    N_t = (up / rho) * dY
    N_p = (up / rho) * (1j * m / st) * Y
    return _sph_to_cart_vec(N_r, N_t, N_p, theta, phi)


# --------------------------------------------------------------------------- #
# Векторизованный эвалуатор M,N (по массиву точек) и трансляция A,B
# --------------------------------------------------------------------------- #
def _dtheta_ynm_arr(n: int, m: int, theta, phi):
    """∂_θ Y_n^m на массиве (θ,φ) (DLMF 14.10.5)."""
    if m < 0:
        am = -m
        return (-1) ** am * np.conj(_dtheta_ynm_arr(n, am, theta, phi))
    ct, st = np.cos(theta), np.sin(theta)
    Pn = sp.lpmv(m, n, ct)
    Pn1 = sp.lpmv(m, n - 1, ct) if (n - 1) >= m else np.zeros_like(ct)
    dP = -((n + m) * Pn1 - n * ct * Pn) / st
    norm = math.sqrt((2 * n + 1) / (4.0 * math.pi)
                     * math.factorial(n - m) / math.factorial(n + m))
    return norm * dP * np.exp(1j * m * phi)


def _sph2cart_arr(Ar, At, Ap, theta, phi):
    """Сферические компоненты (массивы) -> декартовы (Q,3)."""
    ct, st, cp, sp_ = np.cos(theta), np.sin(theta), np.cos(phi), np.sin(phi)
    rhat = np.stack([st * cp, st * sp_, ct], axis=-1)
    that = np.stack([ct * cp, ct * sp_, -st], axis=-1)
    phat = np.stack([-sp_, cp, np.zeros_like(cp)], axis=-1)
    return (Ar[..., None] * rhat + At[..., None] * that + Ap[..., None] * phat)


def mn_grid(n: int, m: int, k: float, pts, kind: str = "reg"):
    """M_{nm}, N_{nm} (декартовы, (Q,3)) на массиве точек pts (Q,3)."""
    r, theta, phi = _to_spherical(pts)
    rho = k * r
    z, up = _zfuncs(n, rho, kind)
    Y = ynm(n, m, theta, phi)
    dY = _dtheta_ynm_arr(n, m, theta, phi)
    imos = (1j * m / np.sin(theta)) * Y
    M = _sph2cart_arr(np.zeros_like(Y), z * imos, -z * dY, theta, phi)
    N = _sph2cart_arr(n * (n + 1) * (z / rho) * Y, (up / rho) * dY, (up / rho) * imos,
                      theta, phi)
    return M, N


def _sphere_quad_cart(R: float, ntheta: int, nphi: int):
    """Точки на сфере радиуса R и веса dΩ (Гаусс–Лежандр × равномерно по φ)."""
    x, wx = np.polynomial.legendre.leggauss(ntheta)
    theta = np.arccos(x)
    phi = 2 * math.pi * np.arange(nphi) / nphi
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    W = np.outer(wx, np.full(nphi, 2 * math.pi / nphi)).ravel()
    TH, PH = TH.ravel(), PH.ravel()
    pts = R * np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], axis=-1)
    return pts, W


def _pick_R(k: float, nmax: int, Rmax: float | None = None) -> float:
    """Радиус сферы проекции, минимизирующий риск нулей j_n (n≤nmax).

    При заданном Rmax (для исходящей трансляции R<|d|) выбор ограничен kR<0.95·k·Rmax.
    """
    hi = 1.4 * nmax + 4.0
    if Rmax is not None:
        hi = min(hi, 0.95 * k * Rmax)
    lo = min(0.6 * nmax + 1.0, 0.5 * hi)
    cand = np.linspace(lo, hi, 300)
    best, bestval = cand[0], -1.0
    for kr in cand:
        v = min(abs(sp.spherical_jn(n, kr)) for n in range(nmax + 1))
        if v > bestval:
            bestval, best = v, kr
    return best / k


def mode_list(nmax: int):
    """Список мод (n,m) для n=1..nmax, m=−n..n (порядок индексации блоков GMM)."""
    return [(n, m) for n in range(1, nmax + 1) for m in range(-n, n + 1)]


def translation_block(d, k: float, nmax: int, source_kind: str = "out"):
    """Полные матрицы векторной трансляции A, B (K×K) и список мод.

    source_kind='reg' — Rg→Rg (j_n в ядре); 'out' — out→Rg (h_n^{(1)}, для GMM:
    рассеянное поле одного тела как падающее для другого). Связь:
        a^M_p = Σ_q (A·c^M_q + B·c^N_q),   a^N_p = Σ_q (B·c^M_q + A·c^N_q).
    """
    d = np.asarray(d, dtype=float)
    dnorm = float(np.linalg.norm(d))
    Rmax = dnorm if source_kind == "out" else None
    R = _pick_R(k, nmax, Rmax)
    nth = 2 * nmax + 6
    pts, W = _sphere_quad_cart(R, nth, nth)
    modes = mode_list(nmax)
    K, Q = len(modes), pts.shape[0]
    TMg = np.empty((K, Q, 3), complex)
    TNg = np.empty((K, Q, 3), complex)
    for i, (n, m) in enumerate(modes):
        TMg[i], TNg[i] = mn_grid(n, m, k, pts, "reg")
    normM = np.einsum("q,iqc->i", W, np.abs(TMg) ** 2)
    normN = np.einsum("q,iqc->i", W, np.abs(TNg) ** 2)
    shifted = pts - d
    SM = np.empty((K, Q, 3), complex)
    for j, (nu, mu) in enumerate(modes):
        SM[j], _ = mn_grid(nu, mu, k, shifted, source_kind)
    A = np.einsum("q,jqc,iqc->ij", W, SM, np.conj(TMg)) / normM[:, None]
    B = np.einsum("q,jqc,iqc->ij", W, SM, np.conj(TNg)) / normN[:, None]
    return A, B, modes


def translation_AB(nu: int, mu: int, d, k: float, nmax: int,
                   kR: float | None = None):
    """Векторные коэффициенты трансляции A^{νμ}_{nm}, B^{νμ}_{nm} (регулярные).

    Тождество:  RgM_{νμ}(r−d) = Σ_{nm} [A·RgM_{nm}(r) + B·RgN_{nm}(r)],
                RgN_{νμ}(r−d) = Σ_{nm} [B·RgM_{nm}(r) + A·RgN_{nm}(r)].
    Считаются спектральной проекцией на сферу радиуса R (точно по квадратуре).
    """
    R = (kR / k) if kR is not None else _pick_R(k, nmax)
    ntheta = 2 * (nmax + nu) + 6
    nphi = 2 * (nmax + nu) + 6
    pts, W = _sphere_quad_cart(R, ntheta, nphi)
    d = np.asarray(d, dtype=float)
    Msrc, _ = mn_grid(nu, mu, k, pts - d, "reg")     # RgM_{νμ}(r−d)
    A, B = {}, {}
    for n in range(1, nmax + 1):
        for m in range(-n, n + 1):
            Mg, Ng = mn_grid(n, m, k, pts, "reg")
            normM = np.sum(W * np.sum(np.abs(Mg) ** 2, axis=1))
            normN = np.sum(W * np.sum(np.abs(Ng) ** 2, axis=1))
            ipM = np.sum(W * np.sum(Msrc * np.conj(Mg), axis=1))
            ipN = np.sum(W * np.sum(Msrc * np.conj(Ng), axis=1))
            A[(n, m)] = ipM / normM
            B[(n, m)] = ipN / normN
    return A, B


# --------------------------------------------------------------------------- #
# Замкнутая (Cruzan) векторная трансляция: аналитический A + точные веса B
# --------------------------------------------------------------------------- #
def _ab_base_term(nu, mu, n, m, p, dr, dth, dph, k, zfun):
    """Общий радиально-угловой множитель T_p (одинаков для A и B)."""
    s = m - mu
    if abs(s) > p:
        return 0.0 + 0.0j, s
    base = (4.0 * math.pi * (1j) ** (n - nu) * (-1) ** m * (-1j) ** p
            * zfun(p, k * dr) * np.conj(ynm(p, s, dth, dph)))
    return base, s


def _A_closed(nu, mu, n, m, dr, dth, dph, k, zfun) -> complex:
    """Замкнутый коэффициент A^{νμ}_{nm} (Wittmann/Cruzan, проверен ~1e-15)."""
    acc = 0.0 + 0.0j
    for p in range(abs(n - nu), n + nu + 1):
        if (n + nu + p) % 2 != 0:
            continue
        base, s = _ab_base_term(nu, mu, n, m, p, dr, dth, dph, k, zfun)
        if base == 0:
            continue
        G = gaunt(n, -m, nu, mu, p, s)
        if G == 0.0:
            continue
        acc += (nu * (nu + 1) + n * (n + 1) - p * (p + 1)) / (2.0 * n * (n + 1)) * base * G
    return acc


def _B_skeleton(nu, mu, n, m, p, dr, dth, dph, k, zfun) -> complex:
    """Скелет B для одного σ=p (нечётная чётность); вес W(n,ν,p) — снаружи."""
    if (n + nu + p) % 2 != 1:
        return 0.0 + 0.0j
    base, s = _ab_base_term(nu, mu, n, m, p, dr, dth, dph, k, zfun)
    if base == 0:
        return 0.0 + 0.0j
    t = wigner_3j(n, nu, p, -m, mu, s) * wigner_3j(n, nu, p - 1, 0, 0, 0)
    if t == 0.0:
        return 0.0 + 0.0j
    return base * t


@lru_cache(maxsize=None)
def _b_weights(nmax: int, k: float) -> tuple:
    """Универсальные веса W(n,ν,σ) для B (не зависят от |d|/направления/типа волны).

    Определяются раз методом наименьших квадратов из проекционных B (на большом разносе,
    где проекция точна) и затем применяются при ЛЮБОМ разносе. Возвращает dict
    {(n,ν): {σ: W}} (кэшируется по (nmax,k)).
    """
    refs = [2.6, 3.4, 4.2]                       # большие |d| (проекция точна)
    dirs = [np.array([0.4, 0.7, 1.6]), np.array([-0.5, 0.9, 1.1]), np.array([0.8, -0.3, 1.4])]
    samples = []
    for dm in refs:
        for dv in dirs:
            d = dm * dv / np.linalg.norm(dv)
            A, B, modes = translation_block(d, k, nmax, "reg")
            samples.append((d, B, modes))
    jn = lambda p, x: sp.spherical_jn(p, x)
    weights: dict = {}
    for n in range(1, nmax + 1):
        for nu in range(1, nmax + 1):
            sig = [p for p in range(abs(n - nu), n + nu + 1) if (n + nu + p) % 2 == 1]
            if not sig:
                continue
            rows, rhs = [], []
            for d, B, modes in samples:
                dr, dth, dph = (float(v) for v in _to_spherical(d))
                for m in range(-n, n + 1):
                    for mu in range(-nu, nu + 1):
                        row = [_B_skeleton(nu, mu, n, m, p, dr, dth, dph, k, jn) for p in sig]
                        if any(abs(x) > 1e-14 for x in row):
                            rows.append(row)
                            rhs.append(B[modes.index((n, m)), modes.index((nu, mu))])
            if not rows:
                continue
            W, *_ = np.linalg.lstsq(np.array(rows, complex), np.array(rhs, complex), rcond=None)
            weights[(n, nu)] = {p: complex(w) for p, w in zip(sig, W)}
    return (weights,)


def translation_block_closed(d, k: float, nmax: int, source_kind: str = "out"):
    """Замкнутая векторная трансляция (Cruzan): A аналитически, B по кэш-весам.

    Без ограничения k·разнос ≳ nmax (в отличие от translation_block). Возвращает (A,B,modes)
    в том же формате (K×K матрицы), пригодные для gmm.
    """
    d = np.asarray(d, dtype=float)
    dr, dth, dph = (float(v) for v in _to_spherical(d))
    zfun = ((lambda p, x: sp.spherical_jn(p, x)) if source_kind == "reg"
            else (lambda p, x: sp.spherical_jn(p, x) + 1j * sp.spherical_yn(p, x)))
    (W,) = _b_weights(nmax, k)
    modes = mode_list(nmax)
    K = len(modes)
    A = np.zeros((K, K), complex)
    B = np.zeros((K, K), complex)
    for i, (n, m) in enumerate(modes):
        for j, (nu, mu) in enumerate(modes):
            A[i, j] = _A_closed(nu, mu, n, m, dr, dth, dph, k, zfun)
            wnn = W.get((n, nu))
            if wnn:
                B[i, j] = sum(wnn[p] * _B_skeleton(nu, mu, n, m, p, dr, dth, dph, k, zfun)
                              for p in wnn)
    return A, B, modes


def translate_regular(coeffs: dict, d, k: float, nmax: int) -> dict:
    """Транслировать разложение регулярных волн {(ν,μ): c} в новый центр (сдвиг d).

    Возвращает {(n,m): c'} c n=0..nmax, где поле, заданное coeffs относительно старого
    центра, выражено относительно нового (смещённого на d).
    """
    out: dict = {}
    for n in range(nmax + 1):
        for m in range(-n, n + 1):
            acc = 0.0 + 0.0j
            for (nu, mu), c in coeffs.items():
                acc += c * alpha_scalar(nu, mu, n, m, d, k)
            out[(n, m)] = acc
    return out
