# SPDX-License-Identifier: MIT
"""Example 9 — Перехват мощности облучателя раскрывом (spill-over efficiency).

Демонстрация модуля :mod:`green_tensor.feed`: интеграл нормированной ДН облучателя
по конусу полуугла theta0 даёт долю мощности первичного источника, попадающую
в раскрыв (линзу), который виден из фазового центра облучателя под этим углом:

    eta_1(theta0) = int_0^theta0 |F|^2 sin d(theta) / int_0^pi |F|^2 sin d(theta)

Связь с геометрией: theta0 = arcsin(a/d) (линза радиуса a, облучатель на
расстоянии d от её центра); theta0 = 90 град — облучатель на поверхности линзы,
канон диссертации: eta_1 элемента Гюйгенса = 7/8 = 0.875.

Запуск:  python3 feed_spillover.py
Выход:   spillover_efficiency.png / .pdf + таблица в stdout
"""
from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))

from green_tensor.feed import (  # noqa: E402
    cone_half_angle,
    dipole_intensity,
    eta1_huygens_closed,
    feed_directivity,
    huygens_intensity,
    spillover_efficiency,
)

FEEDS = (
    ("Элемент Гюйгенса (апертурный)", huygens_intensity),
    ("Элементарный вибратор", dipole_intensity),
)


def main() -> None:
    # --- таблица eta_1(theta0) --------------------------------------------- #
    marks_deg = (5, 10, 20, 30, 45, 60, 90)
    print("Доля мощности облучателя в конусе раскрыва eta_1(theta0):")
    print(f"{'theta0':>8} {'d/a':>7}", *[f"{name.split(' ')[0]:>12}" for name, _ in FEEDS])
    for deg in marks_deg:
        t0 = np.radians(deg)
        d_over_a = 1.0 / np.sin(t0)
        vals = [spillover_efficiency(f2, t0) for _, f2 in FEEDS]
        print(f"{deg:>7}° {d_over_a:>7.2f}", *[f"{v:>12.4f}" for v in vals])

    # --- сверки ------------------------------------------------------------ #
    eta_canon = spillover_efficiency(huygens_intensity, np.pi / 2)
    assert abs(eta_canon - 0.875) < 1e-10          # канон диссертации (полусфера)
    assert abs(eta1_huygens_closed(np.pi / 2) - 0.875) < 1e-12
    D_eg = feed_directivity(huygens_intensity)
    assert abs(D_eg - 3.0) < 1e-6                  # КНД элемента Гюйгенса
    print(f"\nСверка: eta_1_EG(90°) = {eta_canon:.4f} (канон 0.8750), D_EG = {D_eg:.4f}")
    print(f"Геометрия: конус 5° соответствует d/a = {1/np.sin(np.radians(5)):.2f} "
          f"(cone_half_angle(a=1, d={1/np.sin(np.radians(5)):.2f}) = "
          f"{np.degrees(cone_half_angle(1.0, 1/np.sin(np.radians(5)))):.1f}°)")

    # --- график ------------------------------------------------------------ #
    theta0 = np.linspace(0.0, np.pi, 361)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    theta = np.linspace(0.0, np.pi, 721)
    for name, f2 in FEEDS:
        ax1.plot(np.degrees(theta), [f2(t) for t in theta], label=name)
    ax1.set_xlabel(r"$\theta$, град")
    ax1.set_ylabel(r"$|F(\theta)|^2$")
    ax1.set_title("Нормированные ДН облучателей")
    ax1.set_xlim(0, 180)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)

    for name, f2 in FEEDS:
        ax2.plot(np.degrees(theta0), [spillover_efficiency(f2, t) for t in theta0],
                 label=name)
    ax2.axhline(0.875, ls="--", lw=0.8, color="gray")
    ax2.axvline(90, ls="--", lw=0.8, color="gray")
    ax2.plot([90], [0.875], "o", ms=5, color="black",
             label=r"канон: $\eta_1^{ЭГ}(90°)=7/8$")
    t5 = np.radians(5)
    ax2.plot([5], [spillover_efficiency(huygens_intensity, t5)], "s", ms=5,
             color="tab:red", label=r"конус $5°$: $\eta_1=0{,}0057$")
    ax2.set_xlabel(r"полуугол конуса раскрыва $\theta_0$, град")
    ax2.set_ylabel(r"$\eta_1(\theta_0)$")
    ax2.set_title("Перехват мощности раскрывом (spill-over)")
    ax2.set_xlim(0, 180)
    ax2.set_ylim(0, 1.02)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)

    fig.tight_layout()
    here = os.path.dirname(os.path.abspath(__file__))
    for ext in ("png", "pdf"):
        out = os.path.join(here, f"spillover_efficiency.{ext}")
        fig.savefig(out, dpi=200)
        print("Сохранено:", out)


if __name__ == "__main__":
    main()
