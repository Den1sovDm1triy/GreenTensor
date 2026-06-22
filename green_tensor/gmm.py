"""gmm — единый движок сборки сложной геометрии (Generalized Multiparticle Mie).

Потребляет T-матрицы примитивов (любых геометрий, см. scatterer.py) в общем
сферическом базисе ВСВФ и решает самосогласованную систему с трансляционными
теоремами сложения (vswf.translation_block):

    a_p = d_p + Σ_{q≠p} G^{out}_{pq} · T_q a_q ,   c_p = T_p a_p ,

где d_p — разложение внешней падающей волны относительно центра p, а G^{out}_{pq}
переразлагает рассеянное поле тела q как падающее на тело p. Поля и сечения
геометрии-агностичны — зависят только от T-матриц и положений тел.
"""
from __future__ import annotations

import numpy as np

import vswf


def plane_wave_coeffs(k: float, khat, pol, center, nmax: int) -> np.ndarray:
    """Коэффициенты падающей плоской волны E=pol·e^{i k k̂·r} относительно center.

    Возвращает вектор длины 2K (блоки M,N) проекцией на регулярные ВСВФ.
    """
    modes = vswf.mode_list(nmax)
    K = len(modes)
    R = vswf._pick_R(k, nmax)
    pts, W = vswf._sphere_quad_cart(R, 2 * nmax + 6, 2 * nmax + 6)
    kvec = k * np.asarray(khat, dtype=float)
    phase_c = np.exp(1j * np.dot(kvec, np.asarray(center, dtype=float)))
    E = phase_c * np.asarray(pol, dtype=complex)[None, :] * np.exp(1j * (pts @ kvec))[:, None]
    dM = np.empty(K, complex)
    dN = np.empty(K, complex)
    for i, (n, m) in enumerate(modes):
        Mg, Ng = vswf.mn_grid(n, m, k, pts, "reg")
        normM = np.sum(W * np.sum(np.abs(Mg) ** 2, axis=1))
        normN = np.sum(W * np.sum(np.abs(Ng) ** 2, axis=1))
        dM[i] = np.sum(W * np.sum(E * np.conj(Mg), axis=1)) / normM
        dN[i] = np.sum(W * np.sum(E * np.conj(Ng), axis=1)) / normN
    return np.concatenate([dM, dN])


def solve_cluster(scatterers, k: float, khat, pol, nmax: int):
    """Решить кластер: вернуть (a, c, d) формы (P, 2K) — возбуждающие, рассеянные, падающие."""
    P = len(scatterers)
    modes = vswf.mode_list(nmax)
    blk = 2 * len(modes)
    tvecs = [np.asarray(s.t_vector(k, nmax)) for s in scatterers]
    pos = [np.asarray(s.position, dtype=float) for s in scatterers]

    d = np.concatenate([plane_wave_coeffs(k, khat, pol, pos[p], nmax) for p in range(P)])
    Mfull = np.eye(P * blk, dtype=complex)
    for p in range(P):
        for q in range(P):
            if p == q:
                continue
            A, B, _ = vswf.translation_block(pos[p] - pos[q], k, nmax, "out")
            G = np.block([[A, B], [B, A]])          # 2K×2K: (c^M_q,c^N_q) -> a_p
            Mfull[p * blk:(p + 1) * blk, q * blk:(q + 1) * blk] = -G * tvecs[q][None, :]
    a = np.linalg.solve(Mfull, d).reshape(P, blk)
    c = np.array([tvecs[p] * a[p] for p in range(P)])
    return a, c, d.reshape(P, blk)


def _nmax_from_modes(nmodes: int) -> int:
    """Восстановить nmax из числа мод K=nmax(nmax+2)."""
    nmax = 1
    while nmax * (nmax + 2) < nmodes:
        nmax += 1
    return nmax


def scattered_field(scatterers, c, k: float, pts) -> np.ndarray:
    """Полное рассеянное поле кластера E_sca(r) (исходящие ВСВФ всех тел) в точках pts."""
    nmodes = len(c[0]) // 2
    nmax = _nmax_from_modes(nmodes)
    modes = vswf.mode_list(nmax)
    pts = np.asarray(pts, float)
    E = np.zeros((pts.shape[0], 3), dtype=complex)
    for p, s in enumerate(scatterers):
        shifted = pts - np.asarray(s.position, float)
        cp = c[p]
        for i, (n, m) in enumerate(modes):
            Mg, Ng = vswf.mn_grid(n, m, k, shifted, "out")
            E += cp[i] * Mg + cp[nmodes + i] * Ng
    return E


def scattering_cross_section(scatterers, c, k: float, pol, kR_far: float = 300.0) -> float:
    """C_sca кластера интегралом дальнего поля: ∮|E_sca|²R² dΩ / |E_inc|²."""
    nmax = _nmax_from_modes(len(c[0]) // 2)
    R = kR_far / k
    pts, W = vswf._sphere_quad_cart(R, 2 * nmax + 12, 2 * nmax + 12)
    E = scattered_field(scatterers, c, k, pts)
    flux = np.sum(W * np.sum(np.abs(E) ** 2, axis=1)) * R**2
    return float(flux / np.sum(np.abs(np.asarray(pol)) ** 2))


def _forward_amplitude(scatterers, c, k: float, khat, kR: float):
    """Амплитуда рассеяния F вперёд: F = lim R·e^{−ikR}·E_sca(R·k̂)."""
    R = kR / k
    pt = (R * np.asarray(khat, dtype=float)).reshape(1, 3)
    E = scattered_field(scatterers, c, k, pt)[0]
    return R * np.exp(-1j * k * R) * E


def extinction_cross_section(scatterers, c, k: float, khat, pol,
                             kR=(4000.0, 8000.0)) -> float:
    """C_ext кластера (оптическая теорема) с экстраполяцией ближнего поля по 1/(kR)."""
    x1, x2 = 1.0 / kR[0], 1.0 / kR[1]
    F1 = _forward_amplitude(scatterers, c, k, khat, kR[0])
    F2 = _forward_amplitude(scatterers, c, k, khat, kR[1])
    Finf = (F1 * x2 - F2 * x1) / (x2 - x1)          # Ричардсон: убираем O(1/kR)
    pol = np.asarray(pol, dtype=complex)
    return float((4.0 * np.pi / k) * np.imag(np.vdot(pol, Finf))
                 / np.sum(np.abs(pol) ** 2))


def cross_sections(scatterers, c, k: float, khat, pol) -> dict:
    """Сечения кластера: C_sca (дальнее поле), C_ext (опт. теорема), C_abs=C_ext−C_sca."""
    c_sca = scattering_cross_section(scatterers, c, k, pol)
    c_ext = extinction_cross_section(scatterers, c, k, khat, pol)
    return {"c_sca": c_sca, "c_ext": c_ext, "c_abs": c_ext - c_sca}


def isolated_scattered(scatterer, k: float, khat, pol, nmax: int) -> np.ndarray:
    """Рассеянные коэффициенты ОДИНОЧНОГО тела (c = T·d) — для предела расцепления."""
    d = plane_wave_coeffs(k, khat, pol, scatterer.position, nmax)
    return np.asarray(scatterer.t_vector(k, nmax)) * d
