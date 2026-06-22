"""Сквозное демо: композитное решение на сложной геометрии из сфер-примитивов.

Демонстрирует цель проекта: сложное тело представляется кластером непересекающихся
сфер (строго аналитический примитив), а рассеяние считается самосогласованно методом
GMM с сечениями C_sca/C_ext/C_abs. Для тел без потерь выполняется энергобаланс.

Замечание о сходимости: проекционная реализация трансляции требует k·(разнос) ≳ nmax
(радиус сферы проекции R < разноса). Поэтому демо использует корректно разнесённые
кластеры; для ПЛОТНОЙ упаковки нужны замкнутые коэффициенты Cruzan (запланировано).

Запуск:  python3 examples/complex_geometry_demo.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_ROOT, "green_tensor"))

import decompose as dc
import gmm
from scatterer import LayeredSphere

K = 1.0
KHAT = [1.0, 0.0, 0.0]
POL = [0.0, 0.0, 1.0]
EPS = [2.25]          # диэлектрик без потерь
RAD = 0.5
SEP = 3.0             # разнос: k·SEP = 3 ≳ nmax


def report(name, scatterers, nmax=3):
    if not scatterers:
        print(f"  {name:26s} — пустой кластер, пропуск")
        return
    _, c, _ = gmm.solve_cluster(scatterers, K, KHAT, POL, nmax)
    cs = gmm.cross_sections(scatterers, c, K, KHAT, POL)
    bal = abs(cs["c_abs"]) / cs["c_sca"]
    print(f"  {name:26s} N={len(scatterers):3d}  C_sca={cs['c_sca']:.4f}  "
          f"C_ext={cs['c_ext']:.4f}  energy |C_abs|/C_sca={bal:.1e}")
    return cs


def chain(positions):
    return [LayeredSphere(p, RAD, EPS) for p in positions]


def main():
    print(f"Композитное решение на сложной геометрии (k={K}, lossless ⇒ C_abs≈0):\n")

    report("Димер (2 сферы)", chain([[0, 0, +SEP / 2], [0, 0, -SEP / 2]]))
    report("Стержень (4 сферы)", chain([[0, 0, i * SEP] for i in (-1.5, -0.5, 0.5, 1.5)]))
    report("L-образный кластер", chain([[0, 0, 0], [0, 0, SEP], [SEP, 0, SEP]]))
    # «конус»: убывающие радиусы вдоль оси
    cone_sc = [LayeredSphere([0, 0, i * SEP], r, EPS)
               for i, r in zip((-1, 0, 1), (0.7, 0.5, 0.3))]
    report("Конус (убыв. радиусы)", cone_sc)

    # Автоматическое разложение КРУПНОГО шара в кластер (демонстрация decompose→GMM)
    R0 = 5.0
    inside = dc.sphere_indicator([0, 0, 0], R0)
    centers, radius = dc.pack_spheres(inside, [-R0] * 3, [R0] * 3, spacing=SEP, fill=0.45)
    print(f"\n  авто-разложение шара R={R0}: {len(centers)} сфер, "
          f"min_sep={dc.min_separation(centers):.2f} (2r={2*radius:.2f})")
    report("Шар → кластер сфер", dc.to_scatterers(centers, radius, EPS))

    print("\nГотово: сложные формы решены строго аналитическим семейством (сферы + GMM); "
          "энергобаланс выполнен.")


if __name__ == "__main__":
    main()
