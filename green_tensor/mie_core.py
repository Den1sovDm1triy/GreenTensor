"""mie_core — верифицированное ядро рассеяния на многослойной сфере (метод Ми / ТФГ).

Чистая, импортируемая реализация (без matplotlib/input/глобального состояния),
вынесенная из проверенного 01_sphere.py. Содержит ВСЕ исправления, подтверждённые
аналитическим арбитром (tests/analytic_mie.py):

  * внешняя среда = воздух (eps=1), дописывается один раз (бывший баг C4);
  * импедансы по всем гармоникам range(toch) (бывший off-by-one C3);
  * комплексный показатель преломления sigma = sqrt(eps*mu) — поглощение
    входит в аргумент функций Бесселя;
  * фаза eps через atan2 (cmath.phase) — корректно для металлов с Re(eps)<0;
  * экспоненциально масштабированные jve/yve — численная устойчивость для
    больших |eps| (PEC-предел) без переполнения.

Параметры класса:
  polarization : "linear" | "circular"      — тип поляризации падающего поля;
  problem      : "diffraction" | "antenna"  — тип задачи:
        "diffraction" — рассеяние плоской волны (стандартные коэффициенты Ми);
        "antenna"     — возбуждение источником на поверхности сферы
                        («рупор на поверхности», Mn=(Z-1j)/(Z*mH-mHpr)).

Конвенция: e^{-iωt}, исходящая волна ~ h_n^{(1)}. Im(eps) > 0 — поглощение.
Радиусы a[i] нормированы так, что внешний радиус = 1, а k0 — параметр размера k0*R.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import scipy.special as sp

POLARIZATIONS = ("linear", "circular")
PROBLEMS = ("diffraction", "antenna")
ORIENTATIONS = ("vertical", "horizontal")  # для линейной поляризации


# --------------------------------------------------------------------------- #
# Спецфункции Риккати–Бесселя (масштабированные jve/yve, комплексный аргумент)
# --------------------------------------------------------------------------- #
def _psi_chi(n: np.ndarray, z: complex):
    """ψ_n, ψ'_n, χ_n(=Нейман), χ'_n для (комплексного) аргумента z.

    Используются jve/yve (масштаб exp(-|Im z|)) — он сокращается во всех
    отношениях и кросс-произведениях C/S, поэтому при больших |Im z| нет
    переполнения, а результат математически тот же.
    """
    s = np.sqrt(z * np.pi / 2.0)
    a = n / (2 * n + 1)
    b = (n + 1) / (2 * n + 1)
    j_m, j_0, j_p = sp.jve(n - 0.5, z), sp.jve(n + 0.5, z), sp.jve(n + 1.5, z)
    y_m, y_0, y_p = sp.yve(n - 0.5, z), sp.yve(n + 0.5, z), sp.yve(n + 1.5, z)
    psi = s * j_0
    chi = s * y_0
    psi_p = a * s * j_m - b * s * j_p + psi / z
    chi_p = a * s * y_m - b * s * y_p + chi / z
    return psi, psi_p, chi, chi_p


def _xi(n: np.ndarray, z: complex):
    """ξ_n = z·h_n^{(1)} и ξ'_n (для внешней среды; z вещественный → без переполнения)."""
    s = np.sqrt(z * np.pi / 2.0)
    a = n / (2 * n + 1)
    b = (n + 1) / (2 * n + 1)
    h_m, h_0, h_p = sp.hankel1(n - 0.5, z), sp.hankel1(n + 0.5, z), sp.hankel1(n + 1.5, z)
    xi = s * h_0
    xi_p = a * s * h_m - b * s * h_p + xi / z
    return xi, xi_p


@dataclass
class MieSphere:
    """Рассеяние на радиально-слоистой сфере.

    k0   : параметр размера на внешнем радиусе (k0*R, R=1).
    a    : нормированные радиусы границ слоёв (последний = 1, внешний).
    eps  : комплексные диэлектрические проницаемости слоёв.
    miy  : магнитные проницаемости слоёв (по умолчанию 1).
    toch : число гармоник ряда (по умолчанию по правилу Уискомба).
    """

    k0: float
    a: list
    eps: list
    miy: list | None = None
    toch: int | None = None
    polarization: str = "linear"
    problem: str = "diffraction"
    orientation: str = "vertical"   # линейная: "vertical" -> E_theta | "horizontal" -> E_phi
    phi: float = 0.0                # азимут наблюдения для линейной поляризации (E_θ/E_φ)

    _Mn: np.ndarray | None = field(default=None, init=False, repr=False)
    _Nn: np.ndarray | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if self.polarization not in POLARIZATIONS:
            raise ValueError(f"polarization must be one of {POLARIZATIONS}")
        if self.problem not in PROBLEMS:
            raise ValueError(f"problem must be one of {PROBLEMS}")
        if self.orientation not in ORIENTATIONS:
            raise ValueError(f"orientation must be one of {ORIENTATIONS}")
        self.a = [complex(x).real for x in self.a]
        self.eps = [complex(x) for x in self.eps]
        n_layers = len(self.a)
        if len(self.eps) != n_layers:
            raise ValueError("len(eps) must equal len(a)")
        if self.miy is None:
            self.miy = [1.0] * n_layers
        self.miy = [complex(x) for x in self.miy]
        if self.toch is None:
            x = abs(self.k0)
            self.toch = int(math.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))
        self.toch = max(int(self.toch), 1)

    # ------------------------------------------------------------------ ядро
    def coefficients(self):
        """Коэффициенты рассеяния Mn, Nn (кэшируются)."""
        if self._Mn is not None:
            return self._Mn, self._Nn

        nL = len(self.a)
        eps = list(self.eps)
        miy = list(self.miy)
        # внешняя среда (воздух) — ровно один раз, для замыкания рекурсии
        eps_ext = eps + [1.0 + 0j]
        sigma = [np.sqrt(eps[i] * miy[i]) for i in range(nL)]  # компл. показатель

        def kk(j1, j2):  # k[radius j1][layer j2]
            return self.k0 * self.a[j1] * sigma[j2]

        n = np.arange(1, self.toch + 1)

        # внутреннее ядро (слой 0)
        psi0, psi0p, _, _ = _psi_chi(n, kk(0, 0))
        J, Jpr = psi0, psi0p

        Z = np.zeros((self.toch, nL), dtype=complex)
        Y = np.zeros((self.toch, nL), dtype=complex)
        Z[:, 0] = np.sqrt(eps_ext[1] / eps_ext[0]) * (Jpr / J)
        Y[:, 0] = np.sqrt(eps_ext[0] / eps_ext[1]) * (Jpr / J)

        for h in range(1, nL):
            # C,S при j=h-1: A=k[h][h], B=k[h-1][h]
            psiA, psiAp, chiA, chiAp = _psi_chi(n, kk(h, h))
            psiB, psiBp, chiB, chiBp = _psi_chi(n, kk(h - 1, h))
            C = psiA * chiBp - chiA * psiBp
            Cpr = psiAp * chiBp - chiAp * psiBp
            S = chiA * psiB - psiA * chiB
            Spr = chiAp * psiB - psiAp * chiB

            rZ = np.sqrt(eps_ext[h + 1] / eps_ext[h])
            rY = np.sqrt(eps_ext[h] / eps_ext[h + 1])
            Z[:, h] = rZ * (Cpr + Z[:, h - 1] * Spr) / (C + Z[:, h - 1] * S)
            Y[:, h] = rY * (Cpr + Y[:, h - 1] * Spr) / (C + Y[:, h - 1] * S)

        # модифицированные функции на внешней границе (воздух, z = k0)
        mJ, mJpr, _, _ = _psi_chi(n, self.k0)
        mH, mHpr = _xi(n, self.k0)

        Zo, Yo = Z[:, nL - 1], Y[:, nL - 1]
        if self.problem == "diffraction":
            Mn = (Zo * mJ - mJpr) / (Zo * mH - mHpr)
            Nn = (Yo * mJ - mJpr) / (Yo * mH - mHpr)
        else:  # antenna — возбуждение источником на поверхности сферы
            Mn = (Zo - 1j) / (Zo * mH - mHpr)
            Nn = (Yo - 1j) / (Yo * mH - mHpr)

        Mn, Nn = np.conj(Mn), np.conj(Nn)
        self._Mn, self._Nn = Mn, Nn
        return Mn, Nn

    # ----------------------------------------------------------- сечения
    def cross_sections(self) -> dict:
        """Q_sca, Q_ext, Q_abs, Q_back (нормированы на πR²)."""
        Mn, Nn = self.coefficients()
        n = np.arange(1, self.toch + 1)
        x = self.k0
        q_sca = float((2.0 / x**2) * np.sum((2 * n + 1) * (np.abs(Mn) ** 2 + np.abs(Nn) ** 2)))
        q_ext = float((2.0 / x**2) * np.sum((2 * n + 1) * np.real(Mn + Nn)))
        q_back = float((1.0 / x**2) * np.abs(np.sum((2 * n + 1) * ((-1) ** n) * (Mn - Nn))) ** 2)
        return {"q_sca": q_sca, "q_ext": q_ext, "q_abs": q_ext - q_sca, "q_back": q_back}

    # ----------------------------------------------------- угловые функции
    def _angular(self, theta: np.ndarray):
        """π_n, τ_n по углам theta (массивы [toch, len(theta)])."""
        cos_t = np.cos(theta)
        pii = np.zeros((self.toch, len(theta)))
        tay = np.zeros((self.toch, len(theta)))
        for i in range(self.toch):
            m = i + 1
            Lm0 = sp.lpmv(0, m, cos_t)
            Lm1 = sp.lpmv(1, m, cos_t)
            Lm2 = sp.lpmv(2, m, cos_t) if m >= 2 else np.zeros_like(cos_t)
            with np.errstate(divide="ignore", invalid="ignore"):
                base = Lm1 / np.sin(theta)
            for z in range(len(theta)):
                th = theta[z]
                if 0 < th < math.pi:
                    pii[i, z] = base[z]
                elif math.pi < th < 2 * math.pi:
                    pii[i, z] = -base[z]
            tay[i, :] = 0.5 * (Lm2 - m * (m + 1) * Lm0)
        return pii, tay

    def default_theta(self) -> np.ndarray:
        """Сетка углов по умолчанию: 0.01°…360° с шагом 1° (как в 01_sphere)."""
        start = 0.01 * math.pi / 180.0
        step = math.pi / 180.0
        steps = int(((360 - 0.01) * (math.pi / 180.0)) / step) + 1
        return start + step * np.arange(steps)

    # --------------------------------------------------------- диаграмма
    def pattern(self, theta: np.ndarray | None = None) -> dict:
        """Диаграмма рассеяния — ОРИГИНАЛЬНЫЕ формулы GreenTensor.

        linear   -> {"theta", "E_theta", "E_phi", "E", "orientation"}
                    E_θ/E_φ — компоненты поля (формула репо, зависит от phi);
                    orientation: "vertical" -> E_theta, "horizontal" -> E_phi.
        circular -> {"theta", "E_op", "E_kp"}
                    E_op — со-поляризация, E_kp — кросс-поляризация (модуль;
                    формула репо, исправлено real -> abs).
        """
        if theta is None:
            theta = self.default_theta()
        theta = np.asarray(theta, dtype=float)
        Mn, Nn = self.coefficients()
        pii, tay = self._angular(theta)
        n = np.arange(1, self.toch + 1)[:, None]
        w = ((2 * n + 1) / (n * (n + 1))) * ((-1) ** n)   # вес с (-1)^n, как в оригинале
        Mn_c, Nn_c = Mn[:, None], Nn[:, None]

        if self.polarization == "circular":
            # Оригинальная формула репо; real заменён на abs (диаграмма = модуль).
            E_op = np.sum(w * (tay - pii) * (Mn_c + Nn_c), axis=0)
            E_kp = np.sum(w * (tay + pii) * (Mn_c - Nn_c), axis=0)
            return {"theta": theta, "E_op": np.abs(E_op), "E_kp": np.abs(E_kp)}

        # Линейная поляризация — оригинальная формула E_θ/E_φ (зависит от phi).
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
        return {"theta": theta, "E_theta": E_teta, "E_phi": E_phi,
                "E": E, "orientation": self.orientation}
