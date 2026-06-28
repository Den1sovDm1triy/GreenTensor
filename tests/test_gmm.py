# SPDX-License-Identifier: MIT

"""Проверки движка сборки green_tensor/gmm.py (двусферный кластер).

Доказывает корректность всей цепочки связи БЕЗ внешнего бенчмарка:
  [1] предел расцепления: разнос → ∞  ⇒  c_1 → одиночная сфера (c=T·d);
  [2] монотонность: ошибка расцепления убывает с разносом;
  [3] активность связи: на среднем разносе c_1 ЗАМЕТНО ≠ одиночной (связь работает);
  [4] невязка линейной системы решателя ~ 0.

Запуск: python3 tests/test_gmm.py | pytest tests/test_gmm.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _ROOT)

from green_tensor import gmm, vswf  # noqa: E402
from green_tensor.scatterer import LayeredSphere  # noqa: E402

K_BG = 2.0                     # фоновое волновое число
RAD = 0.5                      # радиус сферы (размерный параметр x = k·rad = 1.0)
EPS = [2.25]                   # без потерь
NMAX = 5
KHAT = [1.0, 0.0, 0.0]         # падение вдоль x (поперёк оси кластера z)
POL = [0.0, 0.0, 1.0]          # поляризация вдоль z (⊥ k̂)


def _err_vs_isolated(sep):
    s1 = LayeredSphere([0, 0, +sep / 2], RAD, EPS)
    s2 = LayeredSphere([0, 0, -sep / 2], RAD, EPS)
    _, c, _ = gmm.solve_cluster([s1, s2], K_BG, KHAT, POL, NMAX)
    c_iso = gmm.isolated_scattered(s1, K_BG, KHAT, POL, NMAX)
    return np.max(np.abs(c[0] - c_iso)) / np.max(np.abs(c_iso))


def test_decoupling_limit():
    print("\n[1-2] предел расцепления и монотонность:")
    seps = [4.0 / K_BG, 40.0 / K_BG, 400.0 / K_BG]   # k·s = 4, 40, 400
    errs = [_err_vs_isolated(s) for s in seps]
    for ks, e in zip((4, 40, 400), errs):
        print(f"    k·s={ks:>4}: ||c1 − c_iso||/||c_iso|| = {e:.2e}")
    assert errs[0] > errs[1] > errs[2], "ошибка расцепления должна убывать с разносом"
    assert errs[2] < 5e-3, f"при большом разносе c1 должно сойтись к одиночной: {errs[2]:.1e}"


def test_coupling_is_active():
    print("\n[3] активность связи (средний разнос):")
    e = _err_vs_isolated(8.0 / K_BG)   # k·s = 8 — тела связаны
    print(f"    k·s=8: ||c1 − c_iso||/||c_iso|| = {e:.2e} (должно быть заметным)")
    assert e > 1e-3, "на среднем разносе связь должна заметно менять решение"


def test_solver_residual():
    print("\n[4] невязка линейной системы:")
    s1 = LayeredSphere([0, 0, +1.0], RAD, EPS)
    s2 = LayeredSphere([0, 0, -1.0], RAD, EPS)
    a, c, d = gmm.solve_cluster([s1, s2], K_BG, KHAT, POL, NMAX)
    # проверяем a_p = d_p + Σ_{q≠p} G T a_q напрямую
    A, B, _ = vswf.translation_block_closed(
        s1.position - s2.position, K_BG, NMAX, "out")
    G = np.block([[A, B], [B, A]])
    t2 = s2.t_vector(K_BG, NMAX)
    rhs0 = d[0] + (G * t2[None, :]) @ a[1]
    res = np.max(np.abs(a[0] - rhs0)) / np.max(np.abs(a[0]))
    print(f"    относительная невязка = {res:.1e}")
    assert res < 1e-10


def test_single_sphere_cross_sections():
    """Сечения через gmm (дальнее поле + опт.теорема) == аналитика B&H для одной сферы."""
    import math
    import cmath
    from analytic_mie import q_sca, q_ext
    print("\n[5] сечения кластера (одна сфера) vs B&H:")
    for rad, eps in [(0.5, 2.25), (1.0, 2.25), (0.5, 3.0 + 0.1j)]:
        nmax = max(6, int(K_BG * rad + 4 * (K_BG * rad) ** (1 / 3) + 2))
        s = LayeredSphere([0, 0, 0], rad, [eps])
        c = gmm.isolated_scattered(s, K_BG, KHAT, POL, nmax)
        cs = gmm.cross_sections([s], [c], K_BG, KHAT, POL)
        qsca = cs["c_sca"] / (math.pi * rad**2)
        qext = cs["c_ext"] / (math.pi * rad**2)
        x = K_BG * rad
        d_sca = abs(qsca - q_sca(cmath.sqrt(eps), x, nmax)) / q_sca(cmath.sqrt(eps), x, nmax)
        d_ext = abs(qext - q_ext(cmath.sqrt(eps), x, nmax)) / abs(q_ext(cmath.sqrt(eps), x, nmax))
        print(f"    rad={rad} eps={eps!s:>10}: dQsca={d_sca:.1e} dQext={d_ext:.1e}")
        assert d_sca < 2e-3 and d_ext < 2e-3


def test_energy_conservation_lossless_cluster():
    """Кластер без потерь: C_abs = C_ext − C_sca ≈ 0."""
    print("\n[6] энергобаланс кластера (lossless): C_abs≈0:")
    s1 = LayeredSphere([0, 0, +3.0], 0.5, [2.25])
    s2 = LayeredSphere([0, 0, -3.0], 0.5, [2.25])
    _, c, _ = gmm.solve_cluster([s1, s2], K_BG, KHAT, POL, NMAX)
    cs = gmm.cross_sections([s1, s2], c, K_BG, KHAT, POL)
    rel = abs(cs["c_abs"]) / cs["c_sca"]
    print(f"    C_sca={cs['c_sca']:.4f} C_ext={cs['c_ext']:.4f} |C_abs|/C_sca={rel:.1e}")
    assert rel < 5e-3, f"нарушение энергобаланса кластера: {rel:.1e}"


def test_dense_packing_energy():
    """Плотная упаковка (k·разнос ≪ nmax): замкнутая трансляция сохраняет энергобаланс."""
    import warnings
    print("\n[7] плотная упаковка, энергобаланс (замкнутая Cruzan):")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for sep in (1.2, 1.5):                 # k·sep = 1.2, 1.5 ≪ nmax=3
            pos = [[0, 0, i * sep] for i in (-2, -1, 0, 1, 2)]
            sc = [LayeredSphere(p, 0.5, [2.25]) for p in pos]
            _, c, _ = gmm.solve_cluster(sc, K_BG, KHAT, POL, 3)
            cs = gmm.cross_sections(sc, c, K_BG, KHAT, POL)
            bal = abs(cs["c_abs"]) / cs["c_sca"]
            print(f"    k·sep={K_BG*sep}: C_sca={cs['c_sca']:.4f}  |C_abs|/C_sca={bal:.2e}")
            assert bal < 2e-2, f"энергобаланс плотного кластера нарушен: {bal:.1e}"


def test_fullT_matches_diagonal():
    """Полная (матричная) форма solve_cluster воспроизводит прежнюю поэлементную (диагональную)
    до машинной точности — регрессия Stage 1 (обобщение GMM на недиагональные T-матрицы)."""
    print("\n[8] полная T-матрица == поэлементная диагональная (регрессия):")
    s1 = LayeredSphere([0, 0, +1.0], RAD, EPS)
    s2 = LayeredSphere([0, 0, -1.0], RAD, EPS)
    a_new, c_new, _ = gmm.solve_cluster([s1, s2], K_BG, KHAT, POL, NMAX)
    # эталон: явная ПРЕЖНЯЯ поэлементная сборка (-G * t[None,:], c = t * a)
    pos = [s1.position, s2.position]
    tv = [s1.t_vector(K_BG, NMAX), s2.t_vector(K_BG, NMAX)]
    modes = vswf.mode_list(NMAX)
    blk = 2 * len(modes)
    dref = np.concatenate([gmm.plane_wave_coeffs(K_BG, KHAT, POL, pos[p], NMAX) for p in range(2)])
    M = np.eye(2 * blk, dtype=complex)
    for p in range(2):
        for q in range(2):
            if p == q:
                continue
            A, B, _ = vswf.translation_block_closed(pos[p] - pos[q], K_BG, NMAX, "out")
            G = np.block([[A, B], [B, A]])
            M[p * blk:(p + 1) * blk, q * blk:(q + 1) * blk] = -G * tv[q][None, :]
    aref = np.linalg.solve(M, dref).reshape(2, blk)
    cref = np.array([tv[p] * aref[p] for p in range(2)])
    da = np.max(np.abs(a_new - aref))
    dc = np.max(np.abs(c_new - cref))
    print(f"    max|Δa|={da:.1e}  max|Δc|={dc:.1e}")
    assert da < 1e-12 and dc < 1e-12, "матричная форма обязана совпасть с поэлементной"


class _DiagTSphere:
    """Обёртка-рассеиватель, экспонирующая ПОЛНУЮ T-матрицу t_matrix() = diag(t_vector)
    сферы — проверяет путь недиагональных примитивов (getattr t_matrix)."""

    def __init__(self, sphere):
        self._s = sphere
        self.position = sphere.position

    def t_matrix(self, k, nmax):
        return np.diag(self._s.t_vector(k, nmax))


def test_t_matrix_path_matches_t_vector():
    """Рассеиватель с явной t_matrix == diag(t_vector) даёт тот же кластерный результат,
    что обычная сфера через t_vector — проверка ветки полных T-матриц (для EBCM-примитивов)."""
    print("\n[9] ветка t_matrix (полная) == ветка t_vector (диагональ):")
    s1 = LayeredSphere([0, 0, +1.0], RAD, EPS)
    s2 = LayeredSphere([0, 0, -1.0], RAD, EPS)
    a0, c0, _ = gmm.solve_cluster([s1, s2], K_BG, KHAT, POL, NMAX)
    a1, c1, _ = gmm.solve_cluster([_DiagTSphere(s1), _DiagTSphere(s2)], K_BG, KHAT, POL, NMAX)
    da, dc = np.max(np.abs(a0 - a1)), np.max(np.abs(c0 - c1))
    print(f"    max|Δa|={da:.1e}  max|Δc|={dc:.1e}")
    assert da < 1e-12 and dc < 1e-12, "ветка полной T-матрицы обязана совпасть с диагональной"


def test_axial_incidence_finite():
    """Осевое падение k̂‖ẑ: forward-амплитуда берётся на полюсе θ=0, где раньше
    vswf давал 0/0 ⇒ C_ext/C_abs=NaN. После защиты полюса в mn_grid — конечно;
    для сферы сечения изотропны ⇒ совпадают с поперечным падением."""
    import math
    print("\n[10] осевое падение k̂‖ẑ: C_ext конечно + изотропия сферы:")
    s = LayeredSphere([0, 0, 0], 0.5, [2.25])
    nmax = 6
    c_z = gmm.isolated_scattered(s, K_BG, [0, 0, 1], [1, 0, 0], nmax)   # k̂‖ẑ, pol⊥
    cs_z = gmm.cross_sections([s], [c_z], K_BG, [0, 0, 1], [1, 0, 0])
    print(f"    k̂‖ẑ: C_ext={cs_z['c_ext']:.5f} C_sca={cs_z['c_sca']:.5f} "
          f"|C_abs|/C_sca={abs(cs_z['c_abs'])/cs_z['c_sca']:.1e}")
    assert math.isfinite(cs_z["c_ext"]) and math.isfinite(cs_z["c_abs"]), "осевое падение даёт NaN"
    assert abs(cs_z["c_abs"]) / cs_z["c_sca"] < 5e-3, "lossless: C_abs должно быть ≈0"
    # изотропия сферы: осевое == поперечное падение
    c_x = gmm.isolated_scattered(s, K_BG, [1, 0, 0], [0, 0, 1], nmax)
    cs_x = gmm.cross_sections([s], [c_x], K_BG, [1, 0, 0], [0, 0, 1])
    rel = abs(cs_z["c_ext"] - cs_x["c_ext"]) / cs_x["c_ext"]
    print(f"    изотропия: |C_ext(ẑ)−C_ext(x̂)|/C_ext = {rel:.1e}")
    assert rel < 2e-3, f"сечения сферы должны быть изотропны: {rel:.1e}"


_TESTS = [test_decoupling_limit, test_coupling_is_active, test_solver_residual,
          test_single_sphere_cross_sections, test_energy_conservation_lossless_cluster,
          test_dense_packing_energy, test_fullT_matches_diagonal,
          test_t_matrix_path_matches_t_vector, test_axial_incidence_finite]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ gmm проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
