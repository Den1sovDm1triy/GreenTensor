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


def isolated_scattered(scatterer, k: float, khat, pol, nmax: int) -> np.ndarray:
    """Рассеянные коэффициенты ОДИНОЧНОГО тела (c = T·d) — для предела расцепления."""
    d = plane_wave_coeffs(k, khat, pol, scatterer.position, nmax)
    return np.asarray(scatterer.t_vector(k, nmax)) * d
