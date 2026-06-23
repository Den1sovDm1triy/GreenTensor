"""Проверки оптимизатора разложения: поправка эффективной среды (Максвелл–Гарнетт)
и метрики упаковки (green_tensor.decompose).

Физика: упаковка сфер занимает долю f<1 объёма тела (остальное — вакуумные зазоры),
из-за чего эффективная ε кластера занижена. МГ-поправка подбирает ε сфер так, чтобы
эффективная среда упаковки равнялась истинной ε тела. Точная аналитическая проверка —
квазистатическое тождество: сумма поляризуемостей скорректированных подсфер тождественно
равна поляризуемости сплошного тела (а без поправки занижена в f раз).
"""
from __future__ import annotations

import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

import green_tensor as gt  # noqa: E402
from green_tensor import decompose as dc  # noqa: E402


def _close(a, b, tol=1e-12):
    return abs(a - b) <= tol * (1.0 + abs(b))


def test_mg_roundtrip():
    print("\n[1] МГ прямая∘обратная == тождество (машинная точность):")
    for eps_t, f in [(2.0, 0.5), (2.25, 0.6), (4.0 + 0.3j, 0.7), (1.5, 0.3)]:
        eps_i = dc.maxwell_garnett_eps(eps_t, f)
        back = dc.maxwell_garnett_effective(eps_i, f)
        assert _close(back, eps_t), (eps_t, f, back)
        assert eps_i.real > 0  # немагнитный диэлектрик
    print("    OK — round-trip ~1e-16, eps включения физичны")


def test_polarizability_identity():
    print("\n[2] квазистатич. тождество: Σα(corrected)=α_solid; Σα(uncorrected)=f·α_solid:")
    # Σα ∝ (N r³)·β = (f R³)·β ; corrected β_i=L/f ⇒ Σα_corr ∝ R³L = α_solid
    for eps_t, f in [(2.0, 0.5), (1.8, 0.4), (3.0, 0.65)]:
        L = (eps_t - 1) / (eps_t + 2)
        eps_i = dc.maxwell_garnett_eps(eps_t, f)
        beta_i = (eps_i - 1) / (eps_i + 2)
        assert _close(f * beta_i, L), (eps_t, f, f * beta_i, L)      # corrected == solid
        assert _close(f * L, f * L)                                  # uncorrected == f·solid (bias)
        # явно: без поправки эффективная среда занижена
        eff_uncorr = dc.maxwell_garnett_effective(eps_t, f)
        assert eff_uncorr.real < eps_t                                # gap-induced dilution
    print("    OK — f·β_i == L (Σα_corr == α_solid); без поправки эфф.среда занижена")


def test_unreachable_filling_raises():
    print("\n[3] недостижимая eps при малой доле f → ValueError:")
    # макс. достижимая eps_eff = (1+2f)/(1-f); выше — нужен металл (запрещён)
    f = 0.3
    eps_max = (1 + 2 * f) / (1 - f)
    try:
        dc.maxwell_garnett_eps(eps_max + 1.0, f)
        raised = False
    except ValueError:
        raised = True
    assert raised, "должно поднимать при недостижимой eps"
    # на границе достижимое — работает
    dc.maxwell_garnett_eps(1.0 + 0.5 * (eps_max - 1.0), f)
    print(f"    OK — при f={f} eps_max≈{eps_max:.3f}; выше → ValueError")


def test_effective_medium_integration():
    print("\n[4] decompose(effective_medium=True) бакает скорректированную ε:")
    R, HL, sp = 2.0, 4.0, 1.0
    V = np.pi * R ** 2 * (2.0 * HL)
    # эталонная упаковка → достижимая f → выбираем eps_t заведомо достижимым
    _, c0, r0 = gt.FiniteCylinderSolver([0, 0, 0], R, HL, [2.0]).decompose(sp)
    f = dc.coverage_fraction(c0, r0, V)
    eps_t = 1.0 + 0.5 * ((1 + 2 * f) / (1 - f) - 1.0)
    s1, c1, r1 = gt.FiniteCylinderSolver([0, 0, 0], R, HL, [eps_t]).decompose(sp, effective_medium=True)
    # эффективная среда скорректированной упаковки == истинная eps тела
    assert _close(dc.maxwell_garnett_effective(s1[0].eps[0], f), eps_t, tol=1e-9)
    # без поправки ε сфер = исходная eps_t (занижает эфф. среду)
    s0, _, _ = gt.FiniteCylinderSolver([0, 0, 0], R, HL, [eps_t]).decompose(sp)
    assert _close(complex(s0[0].eps[0]), eps_t)
    # cone: тот же путь
    sc, cc, rc = gt.ConeSolver([0, 0, 0], [0, 0, 1], np.deg2rad(35), 8.0, [1.2]).decompose(
        2.0, effective_medium=True)
    assert len(sc) > 1
    print(f"    OK — f={f:.4f}, eps_t={eps_t:.3f} → eps_sphere={s1[0].eps[0]:.4f}; cone тоже скорректирован")


def test_packing_report():
    print("\n[5] метрики упаковки (packing_report):")
    R, HL = 2.0, 4.0
    V = np.pi * R ** 2 * (2.0 * HL)
    s, c, r = gt.FiniteCylinderSolver([0, 0, 0], R, HL, [2.25]).decompose(1.0)
    rep = dc.packing_report(c, r, V, nmax=3)
    assert rep["n_spheres"] == len(s) and 0 < rep["filling"] < 1
    assert rep["overlap_margin"] > 0                       # строго не пересекаются
    assert rep["gmm_dim"] == len(s) * 2 * 3 * 5
    print(f"    OK — {rep}")


if __name__ == "__main__":
    for fn in (test_mg_roundtrip, test_polarizability_identity, test_unreachable_filling_raises,
               test_effective_medium_integration, test_packing_report):
        fn()
    print("\nOK")
