# SPDX-License-Identifier: MIT

"""nearfield_sphere — полное ближнее поле однородной сферы (метод ТФГ/Ми).

Поле E(r) считается во ВСЕХ точках, включая ВНУТРИ диэлектрика, по диадной функции
Грина слоистой сферы (Li, Kooi, Leong, Yeo, IEEE TMTT 1994): разложение по сферическим
VSWF M_o1n, N_e1n с радиальными функциями j_n (регулярные, внутри) и h_n^{(1)}
(исходящие, снаружи). Внутренние коэффициенты c_n, d_n получены из непрерывности
ТАНГЕНЦИАЛЬНЫХ E и H на границе (поле И производная) — замкнутые формулы Bohren &
Huffman 4.52, эквивалент коэффициентов прохождения метода Li et al. (НЕ деление на
j_n(mka)). Рассеянные a_n, b_n сверяются с каноном 01_sphere по сечению (validate()).

Конвенция: e^{-iωt}, k̂ ∥ z, ê ∥ x (общая ориентация — поворотом точек/поля).
"""
from __future__ import annotations

import math

import numpy as np
import scipy.special as sp


def _riccati(n, z):
    """Риккати-Бессель ψ_n, ψ'_n, ξ_n, ξ'_n. n — массив (Nn,), z — массив (Q,) или скаляр.
    Возвращает массивы формы (Nn, Q) (или (Nn,) для скалярного z)."""
    z = np.asarray(z, dtype=complex)
    scalar = z.ndim == 0
    z = np.atleast_1d(z)
    n = np.asarray(n)[:, None]
    s = np.sqrt(np.pi * z / 2.0)[None, :]
    psi = s * sp.jv(n + 0.5, z[None, :])
    psi_p = s * sp.jv(n - 0.5, z[None, :]) - (n / z[None, :]) * psi
    xi = s * sp.hankel1(n + 0.5, z[None, :])
    xi_p = s * sp.hankel1(n - 0.5, z[None, :]) - (n / z[None, :]) * xi
    if scalar:
        return psi[:, 0], psi_p[:, 0], xi[:, 0], xi_p[:, 0]
    return psi, psi_p, xi, xi_p


def _wiscombe(x):
    return int(np.ceil(abs(x) + 4.0 * abs(x) ** (1.0 / 3.0) + 2.0))


def mie_coeffs(eps, mu, x, nmax):
    """a_n,b_n (рассеяние, исходящие) и c_n,d_n (внутреннее поле, регулярные).

    a_n,b_n — Bohren&Huffman 4.88 (магнитодиэлектрик: m=√(εμ), m̃=√(ε/μ));
    c_n,d_n — B&H 4.52 (непрерывность поля+производной); числитель i·m, знаменатели
    совпадают со знаменателями b_n (для c_n) и a_n (для d_n). При ε=μ=1: a=b=0, c=d=1.
    """
    n = np.arange(1, nmax + 1)
    eps = complex(eps); mu = complex(mu)
    m = np.sqrt(eps * mu); mt = np.sqrt(eps / mu)
    px, ppx, xix, xipx = _riccati(n, x)
    pmx, ppmx, _, _ = _riccati(n, m * x)
    a = (mt * pmx * ppx - px * ppmx) / (mt * pmx * xipx - xix * ppmx)
    b = (pmx * ppx - mt * px * ppmx) / (pmx * xipx - mt * xix * ppmx)
    c = (1j * m) / (pmx * xipx - mt * xix * ppmx)     # TE (M), знам. = знам. b_n
    d = (1j * m) / (mt * pmx * xipx - xix * ppmx)      # TM (N), знам. = знам. a_n
    return n, a, b, c, d


def _pi_tau(nmax, theta):
    """Угловые функции π_n, τ_n (B&H, рекуррентность), форма (nmax, Q)."""
    mu = np.cos(theta)
    Q = mu.size
    pin = np.zeros((nmax, Q)); taun = np.zeros((nmax, Q))
    pin[0] = 1.0
    taun[0] = mu
    if nmax >= 2:
        pin[1] = 3.0 * mu
        taun[1] = 2.0 * mu * pin[1] - 3.0 * pin[0]
    for nn in range(3, nmax + 1):
        i = nn - 1
        pin[i] = ((2 * nn - 1) / (nn - 1)) * mu * pin[i - 1] - (nn / (nn - 1)) * pin[i - 2]
        taun[i] = nn * mu * pin[i] - (nn + 1) * pin[i - 1]
    return pin, taun


def _region_field(A, B, P, Pp, rho, pin, taun, cphi, sphi, sint):
    """Сферические компоненты (Er,Eθ,Eφ) суммы Σ_n [A_n M_o1n + B_n N_e1n].

    z_n=P/ρ; z_n/ρ=P/ρ²; (ρz_n)'/ρ=Pp/ρ.  (всё (nmax,Q); A,B (nmax,1))
    M_o1n: θ=cosφ·π·z, φ=−sinφ·τ·z;  N_e1n: r=cosφ·n(n+1)·sinθ·π·z/ρ,
    θ=cosφ·τ·(ρz)'/ρ, φ=−sinφ·π·(ρz)'/ρ.
    """
    nmax = P.shape[0]
    n = np.arange(1, nmax + 1)[:, None]
    z = P / rho[None, :]
    z_over = P / rho[None, :] ** 2
    rzp = Pp / rho[None, :]
    Er = np.sum(B * (n * (n + 1)) * sint[None, :] * pin * z_over, axis=0) * cphi
    Eth = np.sum(A * pin * z + B * taun * rzp, axis=0) * cphi
    Eph = -np.sum(A * taun * z + B * pin * rzp, axis=0) * sphi
    return Er, Eth, Eph


def field_sphere(a_norm, eps, mu, k, radius, pts, toch=None, khat=(0, 0, 1), pol=(1, 0, 0)):
    """Полное комплексное поле E (Q,3) в точках pts (Q,3): падающее+рассеянное снаружи,
    внутреннее — внутри сферы. Однородная сфера (берётся внешний слой eps[-1] перед
    воздухом? — нет: однородная => eps[0]). Общая ориентация k̂,ê — поворотом."""
    pts = np.asarray(pts, float)
    eps0 = complex(eps[0]); mu0 = complex(mu[0])
    x = k * radius
    nmax = toch or max(_wiscombe(x), 4)
    n, a, b, c, d = mie_coeffs(eps0, mu0, x, nmax)
    En = (1j) ** n * (2 * n + 1) / (n * (n + 1))            # B&H E_n
    m_in = np.sqrt(eps0 * mu0)

    # локальная система: ẑ'=k̂, x̂'=ê⊥, ŷ'=ẑ'×x̂'
    zc = np.asarray(khat, float); zc = zc / (np.linalg.norm(zc) or 1.0)
    xc = np.array(pol, dtype=float)
    xc = xc - np.dot(xc, zc) * zc
    if np.linalg.norm(xc) < 1e-9:
        xc = np.array([1.0, 0, 0]) - zc[0] * zc
    amp = np.linalg.norm(pol) or 1.0
    xc = xc / (np.linalg.norm(xc) or 1.0)
    yc = np.cross(zc, xc)
    R = np.stack([xc, yc, zc], axis=0)                      # глобал->локал (строки = оси)

    loc = pts @ R.T                                         # координаты в локальной системе
    r = np.linalg.norm(loc, axis=1); r = np.where(r < 1e-9, 1e-9, r)
    theta = np.arccos(np.clip(loc[:, 2] / r, -1.0, 1.0))
    phi = np.arctan2(loc[:, 1], loc[:, 0])
    cphi = np.cos(phi); sphi = np.sin(phi); sint = np.sin(theta)
    sint = np.where(np.abs(sint) < 1e-9, 1e-9, sint)
    pin, taun = _pi_tau(nmax, theta)

    Er = np.zeros(pts.shape[0], complex); Eth = np.zeros_like(Er); Eph = np.zeros_like(Er)
    out = r >= radius
    inn = ~out

    if np.any(out):
        ro = r[out]; rho_o = k * ro
        pj, ppj, xh, xph = _riccati(n, rho_o)
        # падающее: A=E_n (M), B=−iE_n (N), регулярные ψ
        Ai = En[:, None]; Bi = (-1j * En)[:, None]
        e1 = _region_field(Ai, Bi, pj, ppj, rho_o, pin[:, out], taun[:, out],
                           cphi[out], sphi[out], sint[out])
        # рассеянное: A=−b_n E_n (M), B=i a_n E_n (N), исходящие ξ
        As = (-b * En)[:, None]; Bs = (1j * a * En)[:, None]
        e2 = _region_field(As, Bs, xh, xph, rho_o, pin[:, out], taun[:, out],
                           cphi[out], sphi[out], sint[out])
        Er[out] = e1[0] + e2[0]; Eth[out] = e1[1] + e2[1]; Eph[out] = e1[2] + e2[2]

    if np.any(inn):
        ri = r[inn]; rho_i = m_in * k * ri
        pj, ppj, _, _ = _riccati(n, rho_i)
        # внутреннее: A=c_n E_n (M), B=−i d_n E_n (N), регулярные ψ(m k r)
        Ac = (c * En)[:, None]; Bc = (-1j * d * En)[:, None]
        e = _region_field(Ac, Bc, pj, ppj, rho_i, pin[:, inn], taun[:, inn],
                          cphi[inn], sphi[inn], sint[inn])
        Er[inn] = e[0]; Eth[inn] = e[1]; Eph[inn] = e[2]

    # (Er,Eθ,Eφ) -> локальные декартовы
    ct = np.cos(theta); st = np.sin(theta)
    ex = Er * st * cphi + Eth * ct * cphi - Eph * sphi
    ey = Er * st * sphi + Eth * ct * sphi + Eph * cphi
    ez = Er * ct - Eth * st
    Eloc = np.stack([ex, ey, ez], axis=1) * amp
    return Eloc @ R                                        # локал->глобал


def _is_metal(eps, mu, k, thickness):
    """Слой «металлический»: Re(eps)<0 (плазмон), очень большой |eps| (PEC-предел),
    или скин-слой меньше толщины (хороший проводник) — поле внутрь не проникает."""
    eps = complex(eps); mu = complex(mu)
    if eps.real < 0.0 or abs(eps) > 1.0e4:
        return True
    n_im = abs(np.sqrt(eps * mu).imag)
    return bool(n_im > 0 and thickness > 0 and 1.0 / (k * n_im) < thickness)


def _assemble(AM, BM, PM, PpM, PN, PpN, rho, pin, taun, cphi, sphi, sint):
    """(Er,Eθ,Eφ) = Σ_n [AM·M_o1n(radial PM) + BN·N_e1n(radial PN)].
    PM/PN — радиальные функции (α ψ+β χ) для TE(M) и TM(N); ρ — аргумент в слое."""
    nmax = PM.shape[0]
    n = np.arange(1, nmax + 1)[:, None]
    zM = PM / rho[None, :]
    zN = PN / rho[None, :]
    zN_over = PN / rho[None, :] ** 2
    rzpN = PpN / rho[None, :]
    Er = np.sum(BM * (n * (n + 1)) * sint[None, :] * pin * zN_over, axis=0) * cphi
    Eth = np.sum(AM * pin * zM + BM * taun * rzpN, axis=0) * cphi
    Eph = -np.sum(AM * taun * zM + BM * pin * rzpN, axis=0) * sphi
    return Er, Eth, Eph


def field_sphere_layered(a_norm, eps, mu, k, radius, pts, toch=None,
                         khat=(0, 0, 1), pol=(1, 0, 0)):
    """Полное поле СЛОИСТОЙ сферы: та же схема, что для однородной, но импедансы Z/Y
    пересчитываются ПО СЛОЯМ каноном 01_sphere. β/α в слое — из Z[l]/Y[l]; амплитуда —
    непрерывностью снаружи внутрь; снаружи — падающее+рассеянное. Металл-слои маскируются
    (поле NaN). Возвращает (E (Q,3) комплекс, mask (Q,) bool — True там, где маска/металл)."""
    pts = np.asarray(pts, float)
    a_norm = [float(x) for x in a_norm]
    epsL = [complex(x) for x in eps]
    muL = [complex(x) for x in mu]
    L = len(a_norm)
    x0 = k * radius
    nmax = toch or max(_wiscombe(x0), 4)
    n = np.arange(1, nmax + 1)

    # параметры слоёв + хост (снаружи)
    sig = [np.sqrt(epsL[l] * muL[l]) for l in range(L)] + [1.0 + 0j]   # относит. показатель
    muH = muL + [1.0 + 0j]
    epH = epsL + [1.0 + 0j]

    def Phi(layer, rnorm, kind):
        """Φ_layer(r): столбцы [ψ,ξ], строки — непрерывные касательные (E,H).
        TE: [U/m; U'/μ];  TM: [U'/m; ε·U/m].  arg = m_layer·x0·rnorm."""
        arg = sig[layer] * x0 * rnorm
        ps, pp, xi, xp = _riccati(n, arg)        # (nmax,)
        m, mu_, ep = sig[layer], muH[layer], epH[layer]
        if kind == "TE":   # непрерывны E∝U/m, H∝U'/μ
            r0 = np.stack([ps / m, xi / m], axis=1)
            r1 = np.stack([pp / mu_, xp / mu_], axis=1)
        else:              # TM (дуально, ψ↔ψ'): E∝U'/m, H∝U/μ
            r0 = np.stack([pp / m, xp / m], axis=1)
            r1 = np.stack([ps / mu_, xi / mu_], axis=1)
        return np.stack([r0, r1], axis=1)        # (nmax,2,2)

    def inv2(M):
        a = M[:, 0, 0]; b = M[:, 0, 1]; c = M[:, 1, 0]; d = M[:, 1, 1]
        det = a * d - b * c
        return np.stack([np.stack([d, -b], axis=1), np.stack([-c, a], axis=1)], axis=1) / det[:, None, None]

    def solve_kind(kind):
        """Трансфер-матрица ядро→хост: коэффициенты [A_l,B_l] слоёв (ψ,ξ) + рассеянный."""
        AB = np.stack([np.ones(nmax, complex), np.zeros(nmax, complex)], axis=1)  # ядро: A=1,B=0
        layers = [AB.copy()]
        for l in range(L):                       # граница a[l] между слоем l и (l+1|хост)
            M = np.einsum('qij,qjk->qik', inv2(Phi(l + 1, a_norm[l], kind)), Phi(l, a_norm[l], kind))
            AB = np.einsum('qij,qj->qi', M, AB)
            layers.append(AB.copy())             # коэффициенты в слое l+1 (последний = хост)
        Ah = AB[:, 0]; Bh = AB[:, 1]
        scale = 1.0 / Ah
        scat = -Bh * scale                       # хост: U = ψ − scat·ξ
        layer_AB = [layers[l] * scale[:, None] for l in range(L)]   # внутр. коэффициенты слоёв
        return scat, layer_AB

    b_n, AB_TE = solve_kind("TE")                # рассеянный TE, коэффициенты слоёв TE
    a_n, AB_TM = solve_kind("TM")               # рассеянный TM, коэффициенты слоёв TM
    px, ppx, xix, xipx = _riccati(n, x0)

    # сборка поля по точкам
    amp = np.linalg.norm(pol) or 1.0
    zc = np.asarray(khat, float); zc = zc / (np.linalg.norm(zc) or 1.0)
    xc = np.array(pol, dtype=float); xc = xc - np.dot(xc, zc) * zc
    if np.linalg.norm(xc) < 1e-9:
        xc = np.array([1.0, 0, 0]) - zc[0] * zc
    xc = xc / (np.linalg.norm(xc) or 1.0); yc = np.cross(zc, xc)
    Rm = np.stack([xc, yc, zc], axis=0)
    loc = pts @ Rm.T
    r = np.linalg.norm(loc, axis=1); r = np.where(r < 1e-9, 1e-9, r)
    rn = r / radius
    theta = np.arccos(np.clip(loc[:, 2] / r, -1, 1)); phi = np.arctan2(loc[:, 1], loc[:, 0])
    cphi = np.cos(phi); sphi = np.sin(phi); sint = np.sin(theta)
    sint = np.where(np.abs(sint) < 1e-9, 1e-9, sint)
    pin, taun = _pi_tau(nmax, theta)
    En = (1j) ** n * (2 * n + 1) / (n * (n + 1))

    Er = np.zeros(pts.shape[0], complex); Eth = np.zeros_like(Er); Eph = np.zeros_like(Er)
    mask = np.zeros(pts.shape[0], bool)
    # какой слой у точки
    bounds = [0.0] + a_norm
    metal_layer = [_is_metal(epsL[l], muL[l], k, (a_norm[l] - bounds[l]) * radius) for l in range(L)]
    # металл экранирует всё внутри: если слой l металл, маскируем r<=a[l]
    metal_to = max([a_norm[l] for l in range(L) if metal_layer[l]] + [0.0])

    # снаружи
    out = rn >= 1.0
    if np.any(out):
        ro = r[out]; rho = k * ro
        pj, ppj, xh, xph = _riccati(n, rho)
        Ai = En[:, None]; Bi = (-1j * En)[:, None]
        e1 = _assemble(Ai, Bi, pj, ppj, pj, ppj, rho, pin[:, out], taun[:, out], cphi[out], sphi[out], sint[out])
        As = (-b_n * En)[:, None]; Bs = (1j * a_n * En)[:, None]
        e2 = _assemble(As, Bs, xh, xph, xh, xph, rho, pin[:, out], taun[:, out], cphi[out], sphi[out], sint[out])
        Er[out] = e1[0] + e2[0]; Eth[out] = e1[1] + e2[1]; Eph[out] = e1[2] + e2[2]

    # по слоям (изнутри): U_l = A_l·ψ(arg) + B_l·ξ(arg), arg = σ_l·x0·r_norm
    for l in range(L):
        sel = (rn >= bounds[l]) & (rn < a_norm[l])
        if not np.any(sel):
            continue
        rho = x0 * rn[sel] * sig[l]
        psi, psip, xi, xip = _riccati(n, rho)
        Ate = AB_TE[l][:, 0][:, None]; Bte = AB_TE[l][:, 1][:, None]
        Atm = AB_TM[l][:, 0][:, None]; Btm = AB_TM[l][:, 1][:, None]
        Ute = Ate * psi + Bte * xi; Upte = Ate * psip + Bte * xip
        Utm = Atm * psi + Btm * xi; Uptm = Atm * psip + Btm * xip
        AM = En[:, None]; BN = (-1j * En)[:, None]
        e = _assemble(AM, BN, Ute, Upte, Utm, Uptm, rho, pin[:, sel], taun[:, sel],
                      cphi[sel], sphi[sel], sint[sel])
        Er[sel] = e[0]; Eth[sel] = e[1]; Eph[sel] = e[2]
        if metal_layer[l]:
            mask[sel] = True
    if metal_to > 0:
        mask |= rn <= metal_to

    ct = np.cos(theta); st = np.sin(theta)
    ex = Er * st * cphi + Eth * ct * cphi - Eph * sphi
    ey = Er * st * sphi + Eth * ct * sphi + Eph * cphi
    ez = Er * ct - Eth * st
    Eloc = np.stack([ex, ey, ez], axis=1) * amp
    return Eloc @ Rm, mask


def validate(eps, mu, k, radius, gt=None):
    """Сверка с каноном 01_sphere по q_sca (если передан gt) + проверка ε=1 → падающее."""
    x = k * radius
    nmax = max(_wiscombe(x), 4)
    n, a, b, c, d = mie_coeffs(eps, mu, x, nmax)
    q_sca = float((2.0 / x ** 2) * np.sum((2 * n + 1) * (np.abs(a) ** 2 + np.abs(b) ** 2)))
    info = {"q_sca": q_sca}
    if gt is not None:
        ref = gt.SphereSolver(radius, [eps], miy=[mu]).cross_sections(k)["q_sca"]
        info["q_sca_core"] = ref
        info["rel_err"] = abs(q_sca - ref) / abs(ref)
    return info
