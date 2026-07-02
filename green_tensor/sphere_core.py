# SPDX-License-Identifier: MIT

"""Adapter for the canonical layered-sphere solver in ``01_sphere.py``.

The mathematical core intentionally stays in ``01_sphere.py``. This module only
loads that file, whose name cannot be imported with normal Python syntax, and
exposes a small importable object used by the public API.
"""
from __future__ import annotations

import importlib.util
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy.special as sp


@lru_cache(maxsize=1)
def _sphere_module():
    path = Path(__file__).with_name("01_sphere.py")
    spec = importlib.util.spec_from_file_location("green_tensor_01_sphere", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load canonical sphere solver from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MieSphere:
    """Importable facade over ``01_sphere.py::RCSCalculator``.

    ``Mn`` and ``Nn`` follow the package convention used by ``01_sphere.py``:
    ``Mn=a_n*`` and ``Nn=b_n*`` for the adopted ``e^{-i omega t}`` convention.
    """

    def __init__(self, k0, a, eps, miy=None, toch=None, **kwargs):
        self.k0 = float(k0)
        self.a = [float(complex(x).real) for x in a]
        self.eps = [complex(x) for x in eps]
        if len(self.a) != len(self.eps):
            raise ValueError("len(a) must equal len(eps)")
        self.miy = [1.0] * len(self.a) if miy is None else [complex(x) for x in miy]
        if len(self.miy) != len(self.a):
            raise ValueError("len(miy) must equal len(a)")
        if toch is None:
            x = abs(self.k0)
            toch = int(math.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
        self.toch = max(int(toch), 1)
        self.polarization = kwargs.get("polarization", "linear")
        self.orientation = kwargs.get("orientation", "vertical")
        self.phi = float(kwargs.get("phi", 0.0))
        self.problem = kwargs.get("problem", "diffraction")
        if self.problem not in ("diffraction", "antenna"):
            raise ValueError("problem must be 'diffraction' or 'antenna'")
        self._coefficients = None

    def _calculator(self):
        calc_cls = _sphere_module().RCSCalculator
        calc = calc_cls(
            k0=self.k0,
            toch=self.toch,
            n=len(self.a),
            phi=1.0,
            a=list(self.a),
            eps=list(self.eps),
            miy=list(self.miy),
            k1=self.k0,
            problem=self.problem,
        )
        calc.run_calculation()
        return calc

    def coefficients(self):
        if self._coefficients is None:
            calc = self._calculator()
            self._coefficients = (np.asarray(calc.Mn, complex), np.asarray(calc.Nn, complex))
        return self._coefficients

    def cross_sections(self) -> dict:
        Mn, Nn = self.coefficients()
        n = np.arange(1, len(Mn) + 1)
        q_sca = float((2.0 / self.k0**2) * np.sum((2 * n + 1) * (np.abs(Mn) ** 2 + np.abs(Nn) ** 2)))
        q_ext = float((2.0 / self.k0**2) * np.sum((2 * n + 1) * np.real(Mn + Nn)))
        q_back = float((1.0 / self.k0**2) * np.abs(np.sum((2 * n + 1) * ((-1) ** n) * (Mn - Nn))) ** 2)
        return {"q_sca": q_sca, "q_ext": q_ext, "q_abs": q_ext - q_sca, "q_back": q_back}

    @staticmethod
    def default_theta() -> np.ndarray:
        return np.arange(0.01, 360.01, 1.0) * math.pi / 180.0

    def _angular(self, theta: np.ndarray):
        cos_t = np.cos(theta)
        pii = np.zeros((self.toch, len(theta)))
        tay = np.zeros((self.toch, len(theta)))
        for i in range(self.toch):
            n = i + 1
            Lm0 = sp.lpmv(0, n, cos_t)
            Lm1 = sp.lpmv(1, n, cos_t)
            Lm2 = sp.lpmv(2, n, cos_t) if n >= 2 else np.zeros_like(cos_t)
            with np.errstate(divide="ignore", invalid="ignore"):
                base = Lm1 / np.sin(theta)
            for z, th in enumerate(theta):
                if 0 < th < math.pi:
                    pii[i, z] = base[z]
                elif math.pi < th < 2 * math.pi:
                    pii[i, z] = -base[z]
            tay[i, :] = 0.5 * (Lm2 - n * (n + 1) * Lm0)
        return pii, tay

    def pattern(self, theta: np.ndarray | None = None) -> dict:
        """Диаграмма рассеяния на сетке углов theta (рад), формулы канона 01_sphere.

        Замечание о нормировочном множителе: для линейной поляризации формула
        решателя содержит norm = (1 - sin²θ·cos²φ)^(-1/2), сингулярный в точках
        θ = 90°/270° при φ = 0 (срез по умолчанию). Сетка по умолчанию
        (default_theta: от 0.01° с шагом 1°) в эти точки не попадает; при задании
        собственной сетки избегайте точных 90°/270° на срезе φ = 0 — там формула
        даёт inf (свойство исходной записи диаграммы канонического решателя).
        """
        if theta is None:
            theta = self.default_theta()
        theta = np.asarray(theta, dtype=float)
        Mn, Nn = self.coefficients()
        pii, tay = self._angular(theta)
        n = np.arange(1, self.toch + 1)[:, None]
        w = ((2 * n + 1) / (n * (n + 1))) * ((-1) ** n)
        Mn_c, Nn_c = Mn[:, None], Nn[:, None]

        if self.polarization == "circular":
            E_op = np.sum(w * (tay - pii) * (Mn_c + Nn_c), axis=0)
            E_kp = np.sum(w * (tay + pii) * (Mn_c - Nn_c), axis=0)
            return {"theta": theta, "E_op": np.abs(E_op), "E_kp": np.abs(E_kp)}

        cphi, sphi = math.cos(self.phi), math.sin(self.phi)
        cos_t = np.cos(theta)
        A = tay * Mn_c - pii * Nn_c
        B = pii * Mn_c - tay * Nn_c
        s_teta = w * (-A * cos_t * cphi ** 2 - B * sphi ** 2)
        s_phi = w * (A * cphi * sphi - B * cos_t ** 2 * sphi * cphi)
        E_teta = np.sum(s_teta, axis=0)
        E_phi = np.sum(s_phi, axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            norm = (1 - (np.sin(theta) * cphi) ** 2) ** (-0.5)
        E_teta = np.abs(norm * E_teta)
        E_phi = np.abs(norm * E_phi)
        E = E_phi if self.orientation == "horizontal" else E_teta
        return {"theta": theta, "E_theta": E_teta, "E_phi": E_phi, "E": E,
                "orientation": self.orientation}
