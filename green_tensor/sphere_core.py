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

from ._angular import pi_tau as _pi_tau


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

    # ---------------------------------------------------------------- #
    # Абсолютные (ненормированные) амплитуды рассеяния и их производные
    # ---------------------------------------------------------------- #
    def amplitudes(self, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        r"""Комплексные амплитудные функции рассеяния S_1(θ), S_2(θ).

        θ — стандартный угол рассеяния (θ=0 — вперёд, θ=π — обратное рассеяние),
        как у Bohren & Huffman.

        .. math::
            S_1(θ) = \sum_n \frac{2n+1}{n(n+1)}\,(M_n \pi_n + N_n \tau_n) \quad(\perp,\ H)
            S_2(θ) = \sum_n \frac{2n+1}{n(n+1)}\,(M_n \tau_n + N_n \pi_n) \quad(\parallel,\ E)

        Здесь ``M_n = a_n^*``, ``N_n = b_n^*`` (конвенция ``e^{-iωt}``); сопряжение
        не влияет на ``|S|²`` и, следовательно, на бистатическую ЭПР и на модуль поля.
        """
        Mn, Nn = self.coefficients()
        theta = np.asarray(theta, dtype=float)
        nmax = len(Mn)
        pin, taun = _pi_tau(nmax, theta)
        n = np.arange(1, nmax + 1)
        w = ((2 * n + 1) / (n * (n + 1)))[:, None]
        Mn_c, Nn_c = Mn[:, None], Nn[:, None]
        S1 = np.sum(w * (Mn_c * pin + Nn_c * taun), axis=0)
        S2 = np.sum(w * (Mn_c * taun + Nn_c * pin), axis=0)
        return S1, S2

    def bistatic_rcs(self, theta: np.ndarray, radius: float) -> dict:
        r"""Абсолютная бистатическая ЭПР σ(θ) в м² по плоскостям наблюдения.

        ``radius`` — внешний физический радиус тела (м). ``self.k0 = x = k·radius``.

        .. math:: \sigma(θ) = \frac{4\pi}{k^2}\,|S(θ)|^2 = \frac{4\pi R^2}{x^2}\,|S(θ)|^2

        Возвращает ``sigma_E`` (плоскость φ=0, параллельная поляризация, S_2) и
        ``sigma_H`` (плоскость φ=π/2, перпендикулярная, S_1). Контроль префактора:
        в точке обратного рассеяния (θ=π) ``sigma_E = sigma_H = π R² · q_back``.
        """
        S1, S2 = self.amplitudes(theta)
        pref = 4.0 * math.pi * float(radius) ** 2 / self.k0**2
        return {"theta": np.asarray(theta, dtype=float),
                "sigma_E": pref * np.abs(S2) ** 2,
                "sigma_H": pref * np.abs(S1) ** 2}

    def far_field(self, theta: np.ndarray, radius: float, e0: float = 1.0,
                  r: float | None = None) -> dict:
        r"""Абсолютная амплитуда рассеянного поля дальней зоны.

        ``k = x/radius`` — физическое волновое число. Величина ``rE = |S|/k·E_0``
        (В) не зависит от расстояния (это ``r·E_scat``) и напрямую сопоставляется с
        HFSS ``rETheta``/``rEPhi``. Если задан ``r`` (м), возвращается также модуль
        поля ``E = |S|/(k r)·E_0`` (В/м) на этом расстоянии.
        """
        S1, S2 = self.amplitudes(theta)
        k = self.k0 / float(radius)
        out = {"theta": np.asarray(theta, dtype=float),
               "rE_E": np.abs(S2) / k * e0,   # φ=0, компонента E_θ
               "rE_H": np.abs(S1) / k * e0}   # φ=π/2, компонента E_φ
        if r is not None:
            out["E_E"] = out["rE_E"] / float(r)
            out["E_H"] = out["rE_H"] / float(r)
        return out

    def field(self, theta: np.ndarray | None = None, e0: complex = 1.0) -> dict:
        r"""Абсолютная напряжённость поля :math:`E_\theta`, :math:`E_\varphi` [В/м]
        с амплитудным множителем источника.

        В теоретической записи диаграммы направленности амплитудный множитель
        поля опускается (диаграмма нормируется). Физическое поле восстанавливается
        умножением безразмерной угловой диаграммы :math:`F(\theta,\varphi)` (метод
        :meth:`pattern`, ур. :math:`E_\theta`/:math:`E_\varphi` главы 2) на
        амплитуду источника (Панченко, ф-ла 3.18):

        .. math::
            \vec E(r,\theta,\varphi) = E_0\,\vec F(\theta,\varphi), \qquad
            E_0 = i\,\frac{I^{e} k_0 l\, Z_0}{4\pi}\,\frac{e^{-i k_0 r}}{r},

        где :math:`E_0` — амплитуда источника, несущая физические единицы и
        радиальный множитель :math:`e^{-ik_0r}/r`. Формула 3.18 отвечает
        элементарному электрическому вибратору (:meth:`dipole_e0`, знаменатель
        :math:`4\pi`, множитель :math:`i`); элемент Гюйгенса задаётся ф-лой 3.14
        (:meth:`huygens_e0`, знаменатель :math:`2\pi` без :math:`i`, по модулю
        вдвое больше).

        .. note::
            Диаграмма направленности линзовой антенны (облучатель на поверхности)
            рассчитывается по соотношению взаимности —
            :func:`green_tensor.nearfield.antenna_pattern_sphere`. Ветвь
            ``problem="antenna"`` отвечает единичному возбуждению поверхности и
            абсолютного уровня облучателя не задаёт.

        Возвращает модули компонент ``E_theta``/``E_phi`` [В/м] на срезе плоскости
        ``self.phi`` и их максимумы ``*_max``.
        """
        pat = self.pattern(theta)
        amp = abs(complex(e0))
        out = {"theta": pat["theta"], "problem": self.problem,
               "polarization": self.polarization, "e0_abs": amp}
        # линейная: E_theta/E_phi/E ; круговая: E_op/E_kp — масштабируем всё, что есть
        for key in ("E_theta", "E_phi", "E", "E_op", "E_kp"):
            if key in pat:
                out[key] = amp * np.asarray(pat[key], dtype=float)
                out[key + "_max"] = float(np.nanmax(out[key]))
        if "orientation" in pat:
            out["orientation"] = pat["orientation"]
        return out

    #: импеданс свободного пространства Z_0 = μ_0 c, Ом
    Z0_FREE_SPACE = 376.730313668

    @staticmethod
    def dipole_e0(current: float, length: float, k0: float, r: float,
                  z0: float = Z0_FREE_SPACE) -> complex:
        r"""Амплитуда поля элементарного вибратора :math:`E_0` (Панченко, ф-ла 3.18).

        .. math:: E_0 = i\,\frac{I^{e} k_0 l\, Z_0}{4\pi}\,\frac{e^{-i k_0 r}}{r}

        ``current`` — ток :math:`I^{e}` [А], ``length`` — длина :math:`l` [м],
        ``k0`` — ФИЗИЧЕСКОЕ волновое число [рад/м] (не электрический размер
        ``self.k0``!), ``r`` — расстояние наблюдения [м], ``z0`` — импеданс
        свободного пространства [Ом]. Возвращает комплексное :math:`E_0` [В/м].
        """
        return 1j * current * k0 * length * z0 / (4.0 * math.pi) * np.exp(-1j * k0 * r) / r

    @staticmethod
    def huygens_e0(current: float, length: float, k0: float, r: float,
                   z0: float = Z0_FREE_SPACE) -> complex:
        r"""Амплитуда поля ЭЛЕМЕНТА ГЮЙГЕНСА :math:`E_0` (Панченко, ф-ла 3.14).

        .. math:: E_0 = \frac{I^{e} k_0 l\, Z_0}{2\pi r}\,e^{-i k_0 r}

        Элемент Гюйгенса — скрещённые тангенциальные электрический
        (:math:`\vec a_\theta I^e`) и магнитный (:math:`\vec a_\varphi I^m/Z_0`) токи;
        его собственная ДН — кардиоида :math:`(1+\cos\theta)/2` (нуль при
        :math:`\theta=\pi`). Именно этот излучатель принят источником в постановке
        задач дифракции и антенной (гл. 2 диссертации).

        Отличается от :meth:`dipole_e0` (ф-ла 3.18, элементарный электрический
        вибратор) знаменателем :math:`2\pi` вместо :math:`4\pi` и отсутствием
        множителя :math:`i`, то есть по модулю вдвое больше.
        """
        return current * k0 * length * z0 / (2.0 * math.pi * r) * np.exp(-1j * k0 * r)

    def circular_field(self, theta: np.ndarray | None, radius: float,
                       e0: complex = 1.0, r: float | None = None) -> dict:
        r"""АБСОЛЮТНЫЕ круговые компоненты рассеянного поля (Панченко, ф-лы 4.14-4.16).

        При падении волны круговой поляризации единичной амплитуды :math:`E_0`

        .. math::
            \vec E^{s} = i\frac{E_0}{\sqrt2}\,\frac{e^{-ik_0r}}{k_0 r}
            \bigl[S_\parallel\vec a_\theta \pm i S_\perp\vec a_\varphi\bigr]e^{i\varphi},

        а проекция на круговые орты :math:`\hat e_{R,L}=(\vec a_\theta\mp i\vec a_\varphi)/\sqrt2`
        даёт

        .. math::
            |E_{R}| = \frac{E_0}{2 k_0 r}\,|S_\parallel-S_\perp|,\qquad
            |E_{L}| = \frac{E_0}{2 k_0 r}\,|S_\parallel+S_\perp|.

        Комбинации совпадают с выдачей :meth:`pattern` при
        ``polarization="circular"``: ``E_op`` :math:`=|S_\parallel-S_\perp|`,
        ``E_kp`` :math:`=|S_\parallel+S_\perp|` (проверено тождественно до
        :math:`10^{-9}`). Множитель :math:`1/2` возникает потому, что величины
        :meth:`pattern` не нормированы на единичную круговую падающую волну.
        Соответствие рукости: ``E_op`` ↔ RHCP, ``E_kp`` ↔ LHCP.

        Возвращает ``E_R``/``E_L`` [В/м] на расстоянии ``r`` (если задано) и
        не зависящие от расстояния ``rE_R``/``rE_L`` [В].
        """
        pat = self.pattern(theta) if self.polarization == "circular" else None
        if pat is None:
            raise ValueError('circular_field требует polarization="circular"')
        k = self.k0 / float(radius)
        half = abs(complex(e0)) / (2.0 * k)
        out = {"theta": pat["theta"],
               "rE_R": half * pat["E_op"], "rE_L": half * pat["E_kp"]}
        if r is not None:
            out["E_R"] = out["rE_R"] / float(r)
            out["E_L"] = out["rE_L"] / float(r)
        return out

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
        r"""Диаграмма рассеяния на сетке углов theta (рад), формулы канона 01_sphere.

        Замечание о нормировочном множителе: для линейной поляризации формула
        решателя содержит norm = (1 - sin²θ·cos²φ)^(-1/2), сингулярный в точках
        θ = 90°/270° при φ = 0 (срез по умолчанию). Сетка по умолчанию
        (default_theta: от 0.01° с шагом 1°) в эти точки не попадает; при задании
        собственной сетки избегайте точных 90°/270° на срезе φ = 0 — там формула
        даёт inf (свойство исходной записи диаграммы канонического решателя).

        Круговая поляризация: ``E_op``/``E_kp`` равны
        :math:`|S_\parallel\mp S_\perp|` и не нормированы на единичную круговую
        падающую волну. Ненормированное поле даёт :meth:`circular_field`
        (множитель :math:`E_0/(2k_0r)`, ф-лы 4.14-4.16 Панченко).
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
