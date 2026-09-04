# SPDX-License-Identifier: MIT

"""Угловые функции π_n, τ_n рассеяния (рекуррентность Бора–Хаффмана).

Чистая (без сингулярного нормировочного множителя) реализация угловых функций,
используемая при сборке амплитудных функций рассеяния S_1, S_2 для АБСОЛЮТНЫХ
диаграмм (:mod:`green_tensor.sphere_core`). Угол θ здесь — стандартный угол
рассеяния (θ=0 — рассеяние вперёд, θ=π — обратное рассеяние), как у Bohren & Huffman.

Independent copy of the standard recurrence kept in the library so the exact core
does not depend on the test-side Mie arbiter (``tests/analytic_mie.py``).
"""
from __future__ import annotations

import numpy as np


def pi_tau(nmax: int, theta: np.ndarray):
    """Массивы π_n, τ_n формы ``[nmax, len(theta)]`` (Bohren & Huffman).

    π_1 = 1;  π_n = ((2n−1)/(n−1))·μ·π_{n−1} − (n/(n−1))·π_{n−2};
    τ_n = n·μ·π_n − (n+1)·π_{n−1};   μ = cos θ.
    """
    mu = np.cos(np.asarray(theta, dtype=float))
    t = mu.shape[0]
    pin = np.zeros((nmax, t))
    taun = np.zeros((nmax, t))
    pin[0] = 1.0
    taun[0] = mu  # 1·μ·π_1 − 2·π_0
    if nmax >= 2:
        pin[1] = 3.0 * mu
        taun[1] = 2.0 * mu * pin[1] - 3.0 * pin[0]
    for n in range(3, nmax + 1):
        i = n - 1
        pin[i] = ((2 * n - 1) / (n - 1)) * mu * pin[i - 1] - (n / (n - 1)) * pin[i - 2]
        taun[i] = n * mu * pin[i] - (n + 1) * pin[i - 1]
    return pin, taun
