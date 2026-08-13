# SPDX-License-Identifier: MIT

"""Строгая валидация многослойного радиально-анизотропного солвера Ми.

Гейты (пороги НЕ ослаблять — при провале чинить реализацию):
  G1  Изотропный предел = многослойное магнитное ядро 01_sphere.py (SphereSolver).
  G2  Один анизотропный слой = anisotropic_mie.py (mie_an_bn_aniso).
  G3  Изотропный предел одиночного слоя = классический Ми (Bohren & Huffman).
  G4  Сохранение энергии: лосслесс -> Q_abs = Q_ext - Q_sca ~ 0.

Плюс демонстрационный расчёт 3-слойной линзы (радиальная перфорация): расхождение
диаграмм рассеяния и Q_back анизотропной модели vs изотропной при k0a=4pi, 8pi.

Запуск:
    cd GreenTensor && PYTHONPATH=. python3 validate_anisotropic_multilayer.py
"""
from __future__ import annotations

import importlib.util
import os
import shutil
import sys
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, "tests"))

from green_tensor.sphere_anisotropic import (  # noqa: E402
    angular_orders,
    bistatic_pattern_anisotropic,
    cross_sections_anisotropic,
    mie_multilayer_anisotropic,
    _riccati_nu,
)
from green_tensor import SphereSolver  # noqa: E402
import analytic_mie as classic  # noqa: E402  (tests/analytic_mie.py — арбитр Bohren&Huffman)


# --------------------------------------------------------------------------- #
# Загрузка эталона одиночной анизотропной сферы (изолированно, вне дерева проекта)
# --------------------------------------------------------------------------- #
def _load_anisotropic_ref():
    """Загрузить anisotropic_mie.py как модуль из изолированной копии.

    При импорте оригинал создаёт рядом каталог out/; чтобы не засорять
    _dissSource, копируем файл во временный каталог и грузим оттуда.
    """
    candidates = [
        os.path.join(os.path.dirname(_HERE), "_dissSource", "codeSource",
                     "python_calc", "anisotropic_sphere", "anisotropic_mie.py"),
    ]
    src = next((p for p in candidates if os.path.exists(p)), None)
    if src is None:
        raise FileNotFoundError("не найден anisotropic_mie.py (эталон одиночного слоя)")
    tmpdir = tempfile.mkdtemp(prefix="aniso_ref_")
    dst = os.path.join(tmpdir, "anisotropic_mie_ref.py")
    shutil.copyfile(src, dst)
    spec = importlib.util.spec_from_file_location("anisotropic_mie_ref", dst)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ANISO_REF = _load_anisotropic_ref()


def _relerr(got, ref):
    got = np.asarray(got); ref = np.asarray(ref)
    denom = np.abs(ref)
    denom = np.where(denom < 1e-300, 1.0, denom)
    return float(np.max(np.abs(got - ref) / denom))


def _abserr(got, ref):
    return float(np.max(np.abs(np.asarray(got) - np.asarray(ref))))


# --------------------------------------------------------------------------- #
# Независимый глобальный солвер (плотная линейная система сшивания полей).
# Совершенно иной алгоритм, чем направленная рекурсия импедансов: собирает
# матрицу 2L×2L из условий непрерывности тангенциальных полей на всех границах
# и решает её напрямую. Условия сшивания (TM): непрерывность Y·u и u' (Y=√(ε_t/μ_t),
# u — Риккати-комбинация слоя, ' — по локальному безразмерному аргументу); TE — дуально.
# Формулировка привязана к валидированному эталону (совпадает с ядром 01_sphere.py
# в изотропном пределе, см. под-проверку в gate_G5).
# --------------------------------------------------------------------------- #
def _global_solve(k0, a_norm, R, eps_r, eps_t, mu_r, mu_t, nmax, mode="TM"):
    a_norm = np.asarray(a_norm, float)
    L = a_norm.size

    def _arr(v):
        v = np.atleast_1d(np.asarray(v, complex))
        return np.full(L, v[0]) if v.size == 1 else v

    eps_r, eps_t, mu_r, mu_t = _arr(eps_r), _arr(eps_t), _arr(mu_r), _arr(mu_t)
    x = k0 * R
    m = np.sqrt(eps_t * mu_t)                       # трансверсальный индекс слоя
    s = x * a_norm                                  # размерные аргументы границ
    Y = np.sqrt(eps_t / mu_t) if mode == "TM" else np.sqrt(mu_t / eps_t)
    pick = 0 if mode == "TM" else 1                 # nu_E или nu_H

    out = np.zeros(nmax, complex)
    for n in range(1, nmax + 1):
        nu = np.array([angular_orders(n, eps_r[l], eps_t[l], mu_r[l], mu_t[l])[pick]
                       for l in range(L)])
        N = 2 * L
        A = np.zeros((N, N), complex)
        rhs = np.zeros(N, complex)
        colA = lambda l: 0 if l == 0 else 1 + 2 * (l - 1)
        colB = lambda l: None if l == 0 else 2 + 2 * (l - 1)
        scat = N - 1
        row = 0
        for l in range(L):
            psi_i, psi_ip, chi_i, chi_ip, _, _ = _riccati_nu(nu[l], m[l] * s[l])
            psi_i, psi_ip = complex(psi_i), complex(psi_ip)
            chi_i, chi_ip = complex(chi_i), complex(chi_ip)
            if l < L - 1:
                lo = l + 1
                psi_o, psi_op, chi_o, chi_op, _, _ = _riccati_nu(nu[lo], m[lo] * s[l])
                psi_o, psi_op = complex(psi_o), complex(psi_op)
                chi_o, chi_op = complex(chi_o), complex(chi_op)
                A[row, colA(l)] += Y[l] * psi_i
                if colB(l) is not None:          # слой 0 регулярен: B_0 отсутствует
                    A[row, colB(l)] += Y[l] * chi_i
                A[row, colA(lo)] -= Y[lo] * psi_o
                A[row, colB(lo)] -= Y[lo] * chi_o
                A[row + 1, colA(l)] += psi_ip
                if colB(l) is not None:
                    A[row + 1, colB(l)] += chi_ip
                A[row + 1, colA(lo)] -= psi_op
                A[row + 1, colB(lo)] -= chi_op
            else:
                pe, pep, _, _, xe, xep = _riccati_nu(float(n), s[l])   # вакуум, порядок n
                pe, pep, xe, xep = complex(pe), complex(pep), complex(xe), complex(xep)
                A[row, colA(l)] += Y[l] * psi_i
                if colB(l) is not None:
                    A[row, colB(l)] += Y[l] * chi_i
                A[row, scat] += xe                       # u_ext = psi_n - scat*xi_n
                rhs[row] = pe
                A[row + 1, colA(l)] += psi_ip
                if colB(l) is not None:
                    A[row + 1, colB(l)] += chi_ip
                A[row + 1, scat] += xep
                rhs[row + 1] = pep
            row += 2
        out[n - 1] = np.linalg.solve(A, rhs)[scat]
    return out


# --------------------------------------------------------------------------- #
# G1 — изотропный предел многослойной линзы vs каноническое ядро (SphereSolver)
# --------------------------------------------------------------------------- #
def gate_G1():
    # 4-слойная градиентная линза (аппрокс. Люнеберга), внешний радиус R=1.
    eps_layers = [1.90, 1.60, 1.30, 1.10]
    a_norm = [0.25, 0.50, 0.75, 1.00]
    R = 1.0
    k0a_list = [0.5, 1.0, 3.0, 5.0, 8.0, 12.0]
    rows = []
    worst = 0.0
    for x in k0a_list:
        k0 = x / R
        nmax = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
        aniso_iso = cross_sections_anisotropic(
            k0, a_norm, R,
            eps_r=eps_layers, eps_t=eps_layers, mu_r=1.0, mu_t=1.0, nmax=nmax)
        core = SphereSolver(R, eps_layers, a_norm=a_norm).cross_sections(k0, toch=nmax)
        d_sca = abs(aniso_iso["Q_sca"] - core["q_sca"]) / abs(core["q_sca"])
        d_back = abs(aniso_iso["Q_back"] - core["q_back"]) / max(abs(core["q_back"]), 1e-300)
        worst = max(worst, d_sca, d_back)
        rows.append((x, aniso_iso["Q_sca"], core["q_sca"], d_sca,
                     aniso_iso["Q_back"], core["q_back"], d_back))
    return worst, rows


# --------------------------------------------------------------------------- #
# G2 — один анизотропный слой vs anisotropic_mie.py (коэффициенты a_n, b_n)
# --------------------------------------------------------------------------- #
def gate_G2():
    # eps_t=1.5, kappa=eps_t/eps_r=0.83 -> eps_r=1.807; mu=1.
    eps_t, eps_r = 1.5, 1.5 / 0.83
    R = 1.0
    cases = [2.0 * np.pi, 4.0 * np.pi, 1.7, 6.0]
    rows = []
    worst = 0.0
    for x in cases:
        nmax = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
        a_me, b_me = mie_multilayer_anisotropic(
            x / R, [1.0], R, eps_r=eps_r, eps_t=eps_t, mu_r=1.0, mu_t=1.0, nmax=nmax)
        a_ref, b_ref = ANISO_REF.mie_an_bn_aniso(
            x, eps_r=eps_r, eps_t=eps_t, mu_r=1.0, mu_t=1.0, n_max=nmax)
        da = _abserr(a_me, a_ref)
        db = _abserr(b_me, b_ref)
        worst = max(worst, da, db)
        rows.append((x, da, db))
    return worst, rows, (eps_r, eps_t)


# --------------------------------------------------------------------------- #
# G3 — изотропный предел одиночного слоя vs классический Ми (Bohren & Huffman)
# --------------------------------------------------------------------------- #
def gate_G3():
    eps = 2.25          # m = 1.5
    m = np.sqrt(eps)
    R = 1.0
    cases = [1.5, 3.0, 2.0 * np.pi]
    rows = []
    worst = 0.0
    for x in cases:
        nmax = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
        a_me, b_me = mie_multilayer_anisotropic(
            x / R, [1.0], R, eps_r=eps, eps_t=eps, mu_r=1.0, mu_t=1.0, nmax=nmax)
        n, a_ref, b_ref = classic.mie_ab(m, x, nmax)
        da = _abserr(a_me, a_ref)
        db = _abserr(b_me, b_ref)
        # плюс сверка сечений
        qs_me = (2.0 / x**2) * np.sum((2 * n + 1) * (np.abs(a_me)**2 + np.abs(b_me)**2))
        qs_ref = classic.q_sca(m, x, nmax)
        dq = abs(qs_me - qs_ref) / abs(qs_ref)
        worst = max(worst, da, db, dq)
        rows.append((x, da, db, dq))
    return worst, rows


# --------------------------------------------------------------------------- #
# G4 — сохранение энергии (лосслесс анизотропный многослойный -> Q_abs ~ 0)
# --------------------------------------------------------------------------- #
def gate_G4():
    eps_r = [1.86, 1.57, 1.28]
    a_norm = [0.53, 0.75, 1.00]
    R = 1.0
    rows = []
    worst = 0.0
    for kappa in [0.70, 0.84, 1.00, 1.25]:
        eps_t = [e * kappa for e in eps_r]      # eps_t/eps_r = kappa
        for x in [1.0, 4.0, 4.0 * np.pi, 8.0 * np.pi]:
            cs = cross_sections_anisotropic(
                x / R, a_norm, R,
                eps_r=eps_r, eps_t=eps_t, mu_r=1.0, mu_t=1.0)
            q_abs = cs["Q_abs"]
            rel = abs(q_abs) / max(abs(cs["Q_sca"]), 1e-300)
            worst = max(worst, rel)
            rows.append((kappa, x, cs["Q_sca"], q_abs, rel))
    return worst, rows


# --------------------------------------------------------------------------- #
# G5 — независимая проверка МНОГОСЛОЙНОЙ АНИЗОТРОПИИ: глобальный солвер vs рекурсия
# --------------------------------------------------------------------------- #
def gate_G5():
    # (a) под-проверка: независимый глобальный солвер сам корректен —
    #     в изотропном пределе совпадает с ядром 01_sphere.py (SphereSolver).
    eps_iso = [1.90, 1.60, 1.30, 1.10]
    a4 = [0.25, 0.50, 0.75, 1.00]
    sub_worst = 0.0
    for x in [1.0, 3.0, 5.0, 8.0]:
        k0 = x
        nmax = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
        aTM = _global_solve(k0, a4, 1.0, eps_iso, eps_iso, 1.0, 1.0, nmax, "TM")
        aTE = _global_solve(k0, a4, 1.0, eps_iso, eps_iso, 1.0, 1.0, nmax, "TE")
        n = np.arange(1, nmax + 1)
        qs = (2.0 / x**2) * np.sum((2 * n + 1) * (np.abs(aTM)**2 + np.abs(aTE)**2))
        core = SphereSolver(1.0, eps_iso, a_norm=a4).cross_sections(k0, toch=nmax)
        sub_worst = max(sub_worst, abs(qs - core["q_sca"]) / abs(core["q_sca"]))

    # (b) ключевая проверка: тот же глобальный солвер vs направленная рекурсия
    #     на МНОГОСЛОЙНОЙ АНИЗОТРОПНОЙ линзе (два aniso-aniso интерфейса).
    eps_r = [1.86, 1.57, 1.28]
    eps_t = [1.56, 1.31, 1.13]
    an = [0.53, 0.75, 1.00]
    rows = []
    worst = 0.0
    for x in [4.0 * np.pi, 8.0 * np.pi, 6.0]:
        k0 = x
        nmax = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
        aTM = _global_solve(k0, an, 1.0, eps_r, eps_t, 1.0, 1.0, nmax, "TM")
        aTE = _global_solve(k0, an, 1.0, eps_r, eps_t, 1.0, 1.0, nmax, "TE")
        a_rec, b_rec = mie_multilayer_anisotropic(k0, an, 1.0, eps_r, eps_t, 1.0, 1.0, nmax)
        da = _abserr(aTM, a_rec)
        db = _abserr(aTE, b_rec)
        worst = max(worst, da, db)
        rows.append((x, da, db))
    return worst, rows, sub_worst


# --------------------------------------------------------------------------- #
# Демонстрация — 3-слойная линза: анизотропная vs изотропная
# --------------------------------------------------------------------------- #
def demonstration():
    eps_r = [1.86, 1.57, 1.28]
    eps_t = [1.56, 1.31, 1.13]
    eps_iso = [1.86, 1.57, 1.28]   # изотропная модель ε = ε_цель
    a_norm = [0.53, 0.75, 1.00]
    R = 1.0
    theta = np.linspace(0.0, np.pi, 721)   # 0..180°, шаг 0.25°
    out = []
    for x in [4.0 * np.pi, 8.0 * np.pi]:
        k0 = x / R
        nmax = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
        P_an = bistatic_pattern_anisotropic(
            k0, a_norm, R, eps_r, eps_t, 1.0, 1.0, theta, nmax=nmax)
        P_iso = bistatic_pattern_anisotropic(
            k0, a_norm, R, eps_iso, eps_iso, 1.0, 1.0, theta, nmax=nmax)
        delta = np.linalg.norm(P_an - P_iso) / np.linalg.norm(P_iso) * 100.0
        cs_an = cross_sections_anisotropic(k0, a_norm, R, eps_r, eps_t, 1.0, 1.0, nmax=nmax)
        cs_iso = cross_sections_anisotropic(k0, a_norm, R, eps_iso, eps_iso, 1.0, 1.0, nmax=nmax)
        dqb = abs(cs_an["Q_back"] - cs_iso["Q_back"]) / abs(cs_iso["Q_back"]) * 100.0
        out.append((x, delta, cs_an["Q_back"], cs_iso["Q_back"], dqb,
                    cs_an["Q_sca"], cs_iso["Q_sca"]))
    return out, (eps_r, eps_t, a_norm)


# --------------------------------------------------------------------------- #
def main() -> int:
    THR = 1e-10
    print("=" * 78)
    print("ВАЛИДАЦИЯ МНОГОСЛОЙНОГО РАДИАЛЬНО-АНИЗОТРОПНОГО СОЛВЕРА МИ")
    print("=" * 78)

    # --- показать, что порядки становятся дробными (санити) ---
    nu_E, nu_H = angular_orders(1, 1.807, 1.5, 1.0, 1.0)
    print(f"\nПроверка дробных порядков (n=1, eps_r=1.807, eps_t=1.5): "
          f"nu_E={nu_E:.6f} (!=1), nu_H={nu_H:.6f} (=1)")
    nu_E2, _ = angular_orders(3, 2.0, 2.0, 1.0, 1.0)
    print(f"Изотропный предел (n=3, eps_r=eps_t): nu_E={nu_E2:.6f} (=3)")

    # --- G1 ---
    g1_worst, g1_rows = gate_G1()
    print("\n" + "-" * 78)
    print("G1  Изотропный предел 4-слойной линзы  vs  ядро SphereSolver (01_sphere.py)")
    print("    eps=[1.90,1.60,1.30,1.10], a=[0.25,0.50,0.75,1.00]")
    print("-" * 78)
    print(f"    {'k0a':>6} {'Q_sca(aniso)':>16} {'Q_sca(core)':>16} {'|d_rel|':>10}"
          f" {'|d_rel|Qb':>10}")
    for x, qsa, qsc, ds, qba, qbc, db in g1_rows:
        print(f"    {x:>6.2f} {qsa:>16.10f} {qsc:>16.10f} {ds:>10.2e} {db:>10.2e}")
    print(f"    worst |d_rel| (Q_sca & Q_back) = {g1_worst:.3e}   "
          f"[{'PASS' if g1_worst < THR else 'FAIL'}]")

    # --- G2 ---
    g2_worst, g2_rows, (er, et) = gate_G2()
    print("\n" + "-" * 78)
    print(f"G2  Один анизотропный слой  vs  anisotropic_mie.py   "
          f"(eps_r={er:.4f}, eps_t={et:.4f}, mu=1)")
    print("-" * 78)
    print(f"    {'k0a':>8} {'max|Δa_n|':>14} {'max|Δb_n|':>14}")
    for x, da, db in g2_rows:
        print(f"    {x:>8.4f} {da:>14.3e} {db:>14.3e}")
    print(f"    worst max|Δ(a_n,b_n)| = {g2_worst:.3e}   "
          f"[{'PASS' if g2_worst < THR else 'FAIL'}]")

    # --- G3 ---
    g3_worst, g3_rows = gate_G3()
    print("\n" + "-" * 78)
    print("G3  Изотропный одиночный слой  vs  классический Ми (Bohren & Huffman), eps=2.25")
    print("-" * 78)
    print(f"    {'k0a':>8} {'max|Δa_n|':>14} {'max|Δb_n|':>14} {'|d_rel|Q_sca':>14}")
    for x, da, db, dq in g3_rows:
        print(f"    {x:>8.4f} {da:>14.3e} {db:>14.3e} {dq:>14.3e}")
    print(f"    worst = {g3_worst:.3e}   [{'PASS' if g3_worst < THR else 'FAIL'}]")

    # --- G4 ---
    g4_worst, g4_rows = gate_G4()
    print("\n" + "-" * 78)
    print("G4  Сохранение энергии (лосслесс анизотропный многослойный): Q_abs = Q_ext - Q_sca")
    print("    3-слойная линза eps_r=[1.86,1.57,1.28], eps_t=kappa*eps_r")
    print("-" * 78)
    print(f"    {'kappa':>6} {'k0a':>8} {'Q_sca':>14} {'Q_abs':>14} {'|Q_abs|/Q_sca':>14}")
    for kappa, x, qs, qa, rel in g4_rows:
        print(f"    {kappa:>6.2f} {x:>8.4f} {qs:>14.8f} {qa:>14.3e} {rel:>14.3e}")
    print(f"    worst |Q_abs|/Q_sca = {g4_worst:.3e}   "
          f"[{'PASS' if g4_worst < THR else 'FAIL'}]")

    # --- G5 ---
    g5_worst, g5_rows, g5_sub = gate_G5()
    print("\n" + "-" * 78)
    print("G5  Многослойная анизотропия: НЕЗАВИСИМЫЙ глобальный солвер (плотная СЛАУ")
    print("    сшивания полей)  vs  направленная рекурсия импедансов — на демо-линзе")
    print("    eps_r=[1.86,1.57,1.28], eps_t=[1.56,1.31,1.13] (два aniso-aniso интерфейса)")
    print("-" * 78)
    print(f"    под-проверка: глобальный солвер vs ядро 01_sphere.py (изотропный предел)"
          f" = {g5_sub:.3e}")
    print(f"    {'k0a':>8} {'max|Δa_n|':>14} {'max|Δb_n|':>14}")
    for x, da, db in g5_rows:
        print(f"    {x:>8.4f} {da:>14.3e} {db:>14.3e}")
    print(f"    worst max|Δ(a_n,b_n)| = {g5_worst:.3e}   "
          f"[{'PASS' if g5_worst < THR else 'FAIL'}]")

    # --- сводка гейтов ---
    gates = {"G1": g1_worst, "G2": g2_worst, "G3": g3_worst, "G4": g4_worst,
             "G5": g5_worst}
    print("\n" + "=" * 78)
    print("ИТОГОВАЯ ТАБЛИЦА ГЕЙТОВ (порог 1e-10)")
    print("=" * 78)
    all_pass = True
    for g, w in gates.items():
        ok = w < THR
        all_pass &= ok
        print(f"    {g}: worst |Δ| = {w:.3e}   ->  {'PASS' if ok else 'FAIL'}")
    print(f"\n    {'ВСЕ ГЕЙТЫ ПРОЙДЕНЫ' if all_pass else 'ЕСТЬ ПРОВАЛЫ — солвер НЕ принят'}")

    if not all_pass:
        print("\nОстановка: демонстрационный расчёт не выполняется до прохождения гейтов.")
        return 1

    # --- демонстрация ---
    demo, (er3, et3, an3) = demonstration()
    print("\n" + "=" * 78)
    print("ДЕМОНСТРАЦИЯ — 3-слойная линза (радиальная перфорация), анизотропная vs изотропная")
    print("=" * 78)
    print(f"    слои: eps_r={er3}, eps_t={et3}  (kappa=eps_t/eps_r)")
    print(f"    kappa по слоям = {[round(t/r, 3) for r, t in zip(er3, et3)]}")
    print(f"    радиусы a_norm={an3}, R=1;  изотропная модель eps=eps_r=[1.86,1.57,1.28]")
    print(f"    δ — относит. расхождение диаграмм рассеяния по норме L2 на θ∈[0,180°]")
    print("-" * 78)
    print(f"    {'k0a':>8} {'δ диаграмм,%':>14} {'Q_back(aniso)':>15} {'Q_back(iso)':>13}"
          f" {'ΔQ_back,%':>11}")
    for x, delta, qba, qbi, dqb, qsa, qsi in demo:
        label = "4π" if abs(x - 4 * np.pi) < 1e-6 else "8π"
        print(f"    {label:>8} {delta:>14.4f} {qba:>15.6f} {qbi:>13.6f} {dqb:>11.4f}")
    print("-" * 78)
    print("    (δ и Q_back посчитаны ЭТИМ провалидированным солвером; iso-ветвь = ε_r=ε_t.)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
