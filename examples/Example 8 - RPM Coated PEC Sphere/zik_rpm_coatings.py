"""
A2: Радиопоглощающие покрытия ZIK на проводящей сфере.

Каноническое применение РПМ — снижение ЭПР металлических объектов.
Параметры взяты verbatim из `_CalculateReferences/ShabuninDifraction/ZIK_BackScat_with_List_ver2.m`:

    f, ГГц = [3,    5.5,    7.5,    10]

    ZIK 51-2190 (4 мм толщиной, диэлектрическое РПМ):
        eps = [3.4-0.03j, 3.4-0.03j, 3.35-0.03j, 3.25-0.03j]
        mu  = [1, 1, 1, 1]

    ZIK 51-2191 (2 мм толщиной, магнетодиэлектрическое РПМ):
        eps = [13.3-1.9j, 12.7-1j, 12.4-0.8j, 12.1-0.7j]
        mu  = [3.1-2.4j, 1.3-2.1j, 0.7-1.6j, 0.44-0.95j]

Метод: PEC-сфера R_pec = 100 мм с однослойным РПМ-покрытием. Вычисляется
эффективность подавления обратного рассеяния (ЭПР) в дБ:

    A(f) = 10 log10[Q_back^bare / Q_back^coated]

Использует:
- mie_pec_sphere() из regression_clean.py для голой PEC-сферы
- mie_pec_with_coating() из multilayer_mie.py для PEC+RPM
"""
from __future__ import annotations
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from green_tensor.multilayer_mie import mie_pec_with_coating
from green_tensor.mie_reference import mie_pec_sphere, cross_sections, wiscombe_nmax

# === ZIK parameters from MATLAB ZIK_BackScat_with_List_ver2.m ==============
F_GHZ = np.array([3.0, 5.5, 7.5, 10.0])
C0 = 2.99792458e8

# ZIK 51-2190 (4 mm, dielectric loss)
EPS_2190 = np.array([3.4-0.03j, 3.4-0.03j, 3.35-0.03j, 3.25-0.03j])
MU_2190 = np.array([1.0, 1.0, 1.0, 1.0], dtype=complex)
T_2190 = 0.004

# ZIK 51-2191 (2 mm, magnetodielectric loss)
EPS_2191 = np.array([13.3-1.9j, 12.7-1j, 12.4-0.8j, 12.1-0.7j])
MU_2191 = np.array([3.1-2.4j, 1.3-2.1j, 0.7-1.6j, 0.44-0.95j])
T_2191 = 0.002

R_PEC = 0.100  # 100 mm


def k0_from_GHz(f_ghz):
    return 2 * np.pi * (f_ghz * 1e9) / C0


def compute_RCS(f_ghz, eps_coat=None, mu_coat=None, t_coat=None):
    """Compute Q_sca, Q_ext, Q_back for PEC sphere R_PEC, optionally with coating."""
    k0 = k0_from_GHz(f_ghz)
    if eps_coat is None:
        # Bare PEC
        x = k0 * R_PEC
        n_max = wiscombe_nmax(x) + 12
        a, b, _ = mie_pec_sphere(x, n_max)
        Qs, Qe, Qb = cross_sections(a, b, x)
        return Qs, Qe, Qb, x
    else:
        # PEC + coating
        R_outer = R_PEC + t_coat
        x = k0 * R_outer
        m_coat = abs(np.sqrt(eps_coat * mu_coat))
        n_max = wiscombe_nmax(m_coat * x) + 12
        a, b = mie_pec_with_coating(k0, R_PEC, [R_outer], [eps_coat], [mu_coat], n_max)
        Qs, Qe, Qb = cross_sections(a, b, x)
        return Qs, Qe, Qb, x


if __name__ == "__main__":
    print("=" * 100)
    print(f"A2: ZIK радиопоглощающие покрытия на PEC-сфере R_PEC = {R_PEC*1000:.0f} мм")
    print("=" * 100)
    print()

    print("Параметры ZIK 51-2190 (4 мм):")
    for i, f in enumerate(F_GHZ):
        print(f"  f={f:5.1f} ГГц: eps={EPS_2190[i]:6.3f}, mu={MU_2190[i]:6.3f}")
    print("\nПараметры ZIK 51-2191 (2 мм):")
    for i, f in enumerate(F_GHZ):
        print(f"  f={f:5.1f} ГГц: eps={EPS_2191[i]:6.3f}, mu={MU_2191[i]:6.3f}")
    print()

    print("=" * 100)
    print(f"{'f, ГГц':>7s} {'k0R':>6s} | {'Q_b bare':>9s} | {'Q_b 2190':>10s} {'A_2190 дБ':>10s} | "
          f"{'Q_b 2191':>10s} {'A_2191 дБ':>10s}")
    print("-" * 100)
    results = {"f_GHz": [], "Qb_bare": [], "Qb_2190": [], "A_2190": [],
               "Qb_2191": [], "A_2191": []}
    for i, f in enumerate(F_GHZ):
        _, _, Qb_bare, x = compute_RCS(f)
        _, _, Qb1, _ = compute_RCS(f, EPS_2190[i], MU_2190[i], T_2190)
        _, _, Qb2, _ = compute_RCS(f, EPS_2191[i], MU_2191[i], T_2191)
        A1 = 10*np.log10(Qb_bare / max(Qb1, 1e-30)) if not np.isnan(Qb1) else float("nan")
        A2 = 10*np.log10(Qb_bare / max(Qb2, 1e-30)) if not np.isnan(Qb2) else float("nan")

        results["f_GHz"].append(f)
        results["Qb_bare"].append(Qb_bare)
        results["Qb_2190"].append(Qb1)
        results["A_2190"].append(A1)
        results["Qb_2191"].append(Qb2)
        results["A_2191"].append(A2)

        Qb2_str = f"{Qb2:10.5e}" if not np.isnan(Qb2) else "       nan"
        A2_str = f"{A2:10.2f}" if not np.isnan(A2) else "       n/a"
        print(f"{f:7.1f} {x:6.2f} | {Qb_bare:9.5f} | {Qb1:10.5f} {A1:10.2f} | "
              f"{Qb2_str} {A2_str}")

    print()
    print("=" * 100)
    print("ИТОГИ (2026-08-13, с высокоточной веткой mpmath в multilayer_mie.py):")
    print("  - ZIK 51-2190 (диэлектрическое, 4 мм): подавление от", end=" ")
    valid_A1 = [a for a in results["A_2190"] if not np.isnan(a)]
    print(f"{min(valid_A1):.1f} до {max(valid_A1):.1f} дБ на 4 частотах")
    valid_A2 = [a for a in results["A_2191"] if not np.isnan(a)]
    if valid_A2:
        print("  - ZIK 51-2191 (магнетодиэлектр., 2 мм): УСИЛЕНИЕ ЭПР на всех 4 частотах:",
              "/".join(f"{-a:.1f}" for a in valid_A2), "дБ")
    print("  - при |Im(m k_0 r)| > 8 расчёт автоматически идёт в arbitrary-precision")
    print("    ветке (mpmath, 60 знаков); double-precision при 3 ГГц даёт ошибку ~1 дБ,")
    print("    выше 5 ГГц — NaN (катастрофическая отмена в лог-производных)")
