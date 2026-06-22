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
sys.path.insert(0, os.path.join(_ROOT, "green_tensor"))

import gmm  # noqa: E402
from scatterer import LayeredSphere  # noqa: E402

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
    A, B, _ = __import__("vswf").translation_block(
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


_TESTS = [test_decoupling_limit, test_coupling_is_active, test_solver_residual,
          test_single_sphere_cross_sections, test_energy_conservation_lossless_cluster]

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
