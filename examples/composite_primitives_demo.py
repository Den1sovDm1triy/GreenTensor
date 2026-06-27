# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Демонстрация: сложное тело как ОБЪЕДИНЕНИЕ аналитических примитивов (CSG) → когерентный GMM.

Собираем «ракетоподобное» тело вдоль оси z из трёх размещённых примитивов, каждый из
которых отдаёт строгую сферическую T-матрицу (сфера — Ми, конус/цилиндр/сфероид — EBCM),
и решаем связанную задачу рассеяния обобщённым методом Ми (GMM):

    конус (нос)  +  конечный цилиндр (корпус)  +  сфероид (хвост)

Печатает сечения кластера (C_sca, C_ext, C_abs) и контроль энергобаланса (для
непоглощающих материалов C_abs≈0). Ось симметрии примитивов ∥ z; положения произвольны.

Запуск:  python3 examples/composite_primitives_demo.py
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from green_tensor import gmm
from green_tensor.scatterer import Cone, FiniteCylinder, Spheroid

K = 2.0                      # волновое число фоновой среды (вакуум)
NMAX = 8
KHAT = (1.0, 0.0, 0.0)      # падение поперёк оси тела
POL = (0.0, 0.0, 1.0)       # поляризация вдоль оси (z)
EPS = 2.25                  # диэлектрик без потерь (стекло-подобный)


def main() -> int:
    body = [
        Cone((0.0, 0.0, 1.6), radius=0.35, height=0.8, eps=EPS),            # нос
        FiniteCylinder((0.0, 0.0, 0.5), radius=0.35, half_length=0.5, eps=EPS),  # корпус
        Spheroid((0.0, 0.0, -0.7), a_eq=0.35, c_ax=0.5, eps=EPS),          # хвост
    ]
    print("Составное тело (CSG-объединение примитивов вдоль z):")
    for p in body:
        print(f"  {type(p).__name__:15s} @ z={p.position[2]:+.2f}  R_bound={p.bounding_radius():.2f}")

    _, c, _ = gmm.solve_cluster(body, K, KHAT, POL, NMAX)
    cs = gmm.cross_sections(body, c, K, KHAT, POL)
    bal = abs(cs["c_abs"]) / cs["c_sca"]
    print(f"\nGMM-решение (k={K}, nmax={NMAX}):")
    print(f"  C_sca = {cs['c_sca']:.5f}")
    print(f"  C_ext = {cs['c_ext']:.5f}")
    print(f"  C_abs = {cs['c_abs']:+.2e}   |C_abs|/C_sca = {bal:.1e}  (lossless ⇒ ≈0)")
    ok = bal < 5e-2
    print("\nЭнергобаланс:", "OK" if ok else "вне допуска (кромки цилиндра/конуса)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
