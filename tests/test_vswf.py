"""Проверки скалярного слоя green_tensor/vswf.py (фундамент GMM).

Самодостаточно (без внешних бенчмарков):
  [1] символы Вигнера 3j — против точных значений и соотношения ортогональности;
  [2] коэффициенты Гаунта — против численной квадратуры ∫YYY dΩ;
  [3] нормировка/ортогональность Y_n^m — против квадратуры;
  [4] скалярная теорема сложения — как ТОЖДЕСТВО ПОЛЕЙ в точках;
  [5] предел нулевого сдвига d→0: α = δ.

Запуск: python3 tests/test_vswf.py | pytest tests/test_vswf.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_ROOT, "green_tensor"))

import vswf  # noqa: E402


def _sphere_quad(nθ=24, nφ=48):
    """Узлы/веса квадратуры по сфере (Гаусс–Лежандр по cosθ, равномерно по φ)."""
    x, wx = np.polynomial.legendre.leggauss(nθ)      # x=cosθ
    theta = np.arccos(x)
    phi = 2 * math.pi * np.arange(nφ) / nφ
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    W = np.outer(wx, np.full(nφ, 2 * math.pi / nφ))  # dΩ = dcosθ dφ
    return TH, PH, W


def test_wigner3j_known():
    print("\n[1] Вигнер 3j:")
    assert abs(vswf.wigner_3j(0, 0, 0, 0, 0, 0) - 1.0) < 1e-14
    assert abs(vswf.wigner_3j(1, 1, 0, 0, 0, 0) - (-1 / math.sqrt(3))) < 1e-14
    assert abs(vswf.wigner_3j(2, 2, 0, 0, 0, 0) - (1 / math.sqrt(5))) < 1e-14
    assert vswf.wigner_3j(1, 1, 1, 0, 0, 0) == 0.0           # сумма j нечётна
    assert vswf.wigner_3j(1, 1, 0, 1, 0, 0) == 0.0           # m1+m2+m3≠0
    # ортогональность (фикс. m3=0): Σ_{m1}(2j3+1)·3j(j1,j2,j3;m1,−m1,0)² = 1
    j1 = j2 = 2
    for j3 in range(abs(j1 - j2), j1 + j2 + 1):
        s = sum((2 * j3 + 1) * vswf.wigner_3j(j1, j2, j3, m1, -m1, 0) ** 2
                for m1 in range(-j1, j1 + 1))
        assert abs(s - 1.0) < 1e-12, f"ортогональность 3j нарушена при j3={j3}: {s}"
    print("    3j: точные значения и ортогональность — ок")


def test_gaunt_vs_quadrature():
    print("\n[2] Гаунт vs квадратура ∫YYY dΩ:")
    TH, PH, W = _sphere_quad()
    cases = [(1, 0, 1, 0, 2, 0), (2, 1, 1, -1, 1, 0), (2, 1, 2, -1, 2, 0),
             (3, 2, 2, -1, 1, -1), (2, 0, 2, 0, 2, 0)]
    for (l1, m1, l2, m2, l3, m3) in cases:
        Y = vswf.ynm(l1, m1, TH, PH) * vswf.ynm(l2, m2, TH, PH) * vswf.ynm(l3, m3, TH, PH)
        num = np.sum(Y * W)
        ana = vswf.gaunt(l1, m1, l2, m2, l3, m3)
        d = abs(num - ana)
        print(f"    G({l1},{m1};{l2},{m2};{l3},{m3}): аналит={ana:+.6e} квадр={num.real:+.6e} d={d:.1e}")
        assert d < 1e-10


def test_ynm_orthonormal():
    print("\n[3] ортонормировка Y_n^m:")
    TH, PH, W = _sphere_quad()
    idx = [(0, 0), (1, 0), (1, 1), (1, -1), (2, 0), (2, 2), (3, -2)]
    for (n, m) in idx:
        for (n2, m2) in idx:
            v = np.sum(vswf.ynm(n, m, TH, PH) * np.conj(vswf.ynm(n2, m2, TH, PH)) * W)
            expect = 1.0 if (n, m) == (n2, m2) else 0.0
            assert abs(v - expect) < 1e-10, f"<Y{n},{m}|Y{n2},{m2}>={v}"
    print("    ⟨Y|Y⟩ = δ — ок")


def test_scalar_addition_theorem():
    """Rg ψ_{νμ}(r−d) = Σ_{nm} α^{νμ}_{nm}(d) Rg ψ_{nm}(r) в точках."""
    print("\n[4] скалярная теорема сложения (тождество полей):")
    k = 1.3
    d = np.array([0.4, 0.9, 1.6])           # |d|≈1.897
    nmax = 14
    pts = [np.array([0.30, -0.20, 0.50]), np.array([-0.40, 0.35, 0.25]),
           np.array([0.10, 0.50, -0.45])]   # |r|<|d| → быстрая сходимость
    # предвычислим α один раз на (ν,μ)
    for (nu, mu) in [(0, 0), (1, 0), (1, 1), (2, -1), (3, 2)]:
        alpha = {(n, m): vswf.alpha_scalar(nu, mu, n, m, d, k)
                 for n in range(nmax + 1) for m in range(-n, n + 1)}
        worst = 0.0
        for r in pts:
            lhs = complex(vswf.rg_scalar(nu, mu, k, r - d))
            rhs = sum(alpha[(n, m)] * complex(vswf.rg_scalar(n, m, k, r))
                      for n in range(nmax + 1) for m in range(-n, n + 1))
            worst = max(worst, abs(lhs - rhs))
        print(f"    (ν,μ)=({nu},{mu}): max|LHS−RHS|={worst:.1e}")
        assert worst < 1e-8, f"теорема сложения не сходится для (ν,μ)=({nu},{mu})"


def test_zero_translation_identity():
    print("\n[5] предел d=0: α = δ:")
    k = 1.7
    z = np.zeros(3)
    assert abs(vswf.alpha_scalar(2, 1, 2, 1, z, k) - 1.0) < 1e-12
    assert abs(vswf.alpha_scalar(2, 1, 3, 1, z, k)) < 1e-12
    assert abs(vswf.alpha_scalar(2, 1, 2, 0, z, k)) < 1e-12
    assert abs(vswf.alpha_scalar(0, 0, 0, 0, z, k) - 1.0) < 1e-12
    print("    α(d=0) = δ — ок")


def _num_curl_div(F, p, h=1e-4):
    """Численные ∇×F и ∇·F (центральные разности) для вектор-функции F(point)->(3,)."""
    e = np.eye(3)
    d = [(np.asarray(F(p + h * e[i])) - np.asarray(F(p - h * e[i]))) / (2 * h) for i in range(3)]
    # d[i] = ∂F/∂x_i ; компоненты d[i][j] = ∂F_j/∂x_i
    curl = np.array([d[1][2] - d[2][1], d[2][0] - d[0][2], d[0][1] - d[1][0]])
    div = d[0][0] + d[1][1] + d[2][2]
    return curl, div


def test_vsw_curl_relations():
    """∇×M=kN, ∇×N=kM, ∇·M=∇·N=0, M·r̂=0 (численная проверка эвалуатора)."""
    print("\n[6] векторные M,N: соотношения ротора/дивергенции:")
    k = 1.2
    pts = [np.array([0.50, 0.40, 0.60]), np.array([-0.45, 0.35, 0.55])]
    for kind in ("reg", "out"):
        for (n, m) in [(1, 0), (1, 1), (2, -1), (3, 2)]:
            for p in pts:
                M = lambda q: vswf.vsw_M(n, m, k, q, kind)
                N = lambda q: vswf.vsw_N(n, m, k, q, kind)
                curlM, divM = _num_curl_div(M, p)
                curlN, divN = _num_curl_div(N, p)
                kN = k * np.asarray(N(p))
                kM = k * np.asarray(M(p))
                e_curlM = np.max(np.abs(curlM - kN)) / np.max(np.abs(kN))
                e_curlN = np.max(np.abs(curlN - kM)) / np.max(np.abs(kM))
                # M·r̂ = 0 (нет радиальной компоненты)
                rhat = p / np.linalg.norm(p)
                m_rad = abs(np.dot(np.asarray(M(p)), rhat))
                assert e_curlM < 1e-5, f"∇×M≠kN ({kind} n={n} m={m}): {e_curlM:.1e}"
                assert e_curlN < 1e-5, f"∇×N≠kM ({kind} n={n} m={m}): {e_curlN:.1e}"
                assert abs(divM) < 1e-4 and abs(divN) < 1e-4, "∇·M или ∇·N ≠ 0"
                assert m_rad < 1e-12, f"M·r̂≠0: {m_rad:.1e}"
        print(f"    {kind}: ∇×M=kN, ∇×N=kM, ∇·=0, M⊥r̂ — ок")


def test_vector_translation():
    """Векторная теорема сложения: RgM_{νμ}(r−d)=Σ[A·RgM+B·RgN], и N-источник."""
    print("\n[7] векторная трансляция A,B (тождество полей, оба источника):")
    k = 1.1
    d = np.array([0.3, 0.7, 1.5])           # |d|≈1.69
    nmax = 10
    pts = [np.array([0.25, -0.15, 0.40]), np.array([-0.30, 0.30, 0.20])]  # |r|<|d|
    for (nu, mu) in [(1, 0), (1, 1), (2, -1)]:
        A, B = vswf.translation_AB(nu, mu, d, k, nmax)
        worstM = worstN = 0.0
        for r in pts:
            lhsM = np.asarray(vswf.vsw_M(nu, mu, k, r - d, "reg"))
            lhsN = np.asarray(vswf.vsw_N(nu, mu, k, r - d, "reg"))
            rhsM = np.zeros(3, complex)
            rhsN = np.zeros(3, complex)
            for n in range(1, nmax + 1):
                for m in range(-n, n + 1):
                    Mr = np.asarray(vswf.vsw_M(n, m, k, r, "reg"))
                    Nr = np.asarray(vswf.vsw_N(n, m, k, r, "reg"))
                    a, b = A[(n, m)], B[(n, m)]
                    rhsM += a * Mr + b * Nr
                    rhsN += b * Mr + a * Nr
            worstM = max(worstM, np.max(np.abs(lhsM - rhsM)))
            worstN = max(worstN, np.max(np.abs(lhsN - rhsN)))
        print(f"    (ν,μ)=({nu},{mu}): max|ΔM|={worstM:.1e}  max|ΔN|={worstN:.1e}")
        assert worstM < 1e-5 and worstN < 1e-5, f"векторная трансляция не сходится ({nu},{mu})"


_TESTS = [test_wigner3j_known, test_gaunt_vs_quadrature, test_ynm_orthonormal,
          test_scalar_addition_theorem, test_zero_translation_identity,
          test_vsw_curl_relations, test_vector_translation]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ vswf проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
