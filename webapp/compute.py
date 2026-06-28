# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Расчётный мост веб-приложения GreenTensor.

Принимает параметры слоёв (JSON), считает через каноническое ядро
green_tensor/sphere_core.py (фасад точного 01_sphere.py) и возвращает серии данных
для диаграмм: линейная/круговая поляризация, обе задачи (дифракция/антенна),
сечения, свип по k0, коэффициенты Mn/Nn.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _ROOT)

from green_tensor.sphere_core import MieSphere  # noqa: E402

DB_FLOOR = -60.0
MAX_TOCH = 150
MAX_SWEEP_POINTS = 400
MAX_LAYERS = 24
# Сетка углов: без точных 0/180/360 (особые точки угловых функций pi_n)
THETA_DEG = np.arange(0.25, 360.0, 0.5)
THETA_RAD = np.deg2rad(THETA_DEG)


class ComputeError(ValueError):
    """Ошибка валидации входных параметров (сообщение показывается в UI)."""


def _finite_or_none(arr):
    """np-массив -> список (не-конечные значения -> None для JSON)."""
    out = []
    for v in np.asarray(arr, dtype=float).ravel():
        out.append(float(v) if np.isfinite(v) else None)
    return out


def _db_pair(e1: np.ndarray, e2: np.ndarray):
    """Совместная dB-нормировка пары кривых к ОБЩЕМУ максимуму (видны отн. уровни)."""
    stack = np.concatenate([np.abs(e1), np.abs(e2)])
    stack = stack[np.isfinite(stack)]
    emax = float(stack.max()) if stack.size else 0.0
    if emax <= 0:
        z = np.full_like(np.asarray(e1, dtype=float), DB_FLOOR)
        return z, z.copy(), emax
    with np.errstate(divide="ignore", invalid="ignore"):
        d1 = 20.0 * np.log10(np.abs(e1) / emax)
        d2 = 20.0 * np.log10(np.abs(e2) / emax)
    return np.maximum(d1, DB_FLOOR), np.maximum(d2, DB_FLOOR), emax


def _parse_layers(payload: dict):
    raw = payload.get("layers")
    if not isinstance(raw, list) or not raw:
        raise ComputeError("Не заданы слои.")
    if len(raw) > MAX_LAYERS:
        raise ComputeError(f"Слишком много слоёв (макс. {MAX_LAYERS}).")
    layers = []
    for i, L in enumerate(raw, 1):
        try:
            r = float(L["r"])
            eps = complex(float(L.get("eps_re", 1.0)), float(L.get("eps_im", 0.0)))
            mu = complex(float(L.get("mu_re", 1.0)), float(L.get("mu_im", 0.0)))
        except (KeyError, TypeError, ValueError):
            raise ComputeError(f"Слой {i}: некорректные значения.")
        if not (r > 0 and math.isfinite(r)):
            raise ComputeError(f"Слой {i}: радиус должен быть > 0.")
        layers.append((r, eps, mu))
    radii = [r for r, _, _ in layers]
    if any(b <= a for a, b in zip(radii, radii[1:])):
        raise ComputeError("Радиусы слоёв должны строго возрастать.")
    r_out = radii[-1]
    # Нормировка к внешнему радиусу + канонический внешний слой-воздух (как в эталоне)
    a = [r / r_out for r, _, _ in layers] + [1.0]
    eps = [e for _, e, _ in layers] + [1.0 + 0j]
    miy = [m for _, _, m in layers] + [1.0 + 0j]
    return a, eps, miy


def _coeff_block(sphere: MieSphere):
    Mn, Nn = sphere.coefficients()
    n = list(range(1, len(Mn) + 1))
    return {"n": n, "Mn_abs": _finite_or_none(np.abs(Mn)), "Nn_abs": _finite_or_none(np.abs(Nn))}


def _problem_block(k0, a, eps, miy, toch, phi_rad, problem):
    base = dict(k0=k0, a=list(a), eps=list(eps), miy=list(miy), toch=toch, problem=problem)
    # Главные плоскости — как в эталонном использовании (01_sphere __main__,
    # сравнения с HFSS): E-плоскость phi=0 и H-плоскость phi=90, обе — E_theta.
    # В оригинальной формуле E_phi ~ sin(phi)cos(phi) и тождественно 0 в главных
    # плоскостях, поэтому пользовательский phi — дополнительный срез.
    lin_e = MieSphere(**base, polarization="linear", phi=0.0)
    lin_h = MieSphere(**base, polarization="linear", phi=math.pi / 2)
    lin_c = MieSphere(**base, polarization="linear", phi=phi_rad)
    cir = MieSphere(**base, polarization="circular")

    pe = lin_e.pattern(THETA_RAD)["E_theta"]
    ph = lin_h.pattern(THETA_RAD)["E_theta"]
    pc = lin_c.pattern(THETA_RAD)
    p_cir = cir.pattern(THETA_RAD)

    db_e, db_h, _ = _db_pair(pe, ph)
    db_ct, db_cp, _ = _db_pair(pc["E_theta"], pc["E_phi"])
    db_op, db_kp, _ = _db_pair(p_cir["E_op"], p_cir["E_kp"])

    block = {
        "linear": {
            "E_plane": {"E": _finite_or_none(pe), "dB": _finite_or_none(db_e)},
            "H_plane": {"E": _finite_or_none(ph), "dB": _finite_or_none(db_h)},
            "custom": {
                "phi_deg": math.degrees(phi_rad),
                "E_theta": _finite_or_none(pc["E_theta"]),
                "E_phi": _finite_or_none(pc["E_phi"]),
                "dB_theta": _finite_or_none(db_ct),
                "dB_phi": _finite_or_none(db_cp),
            },
        },
        "circular": {
            "E_op": _finite_or_none(p_cir["E_op"]),
            "E_kp": _finite_or_none(p_cir["E_kp"]),
            "dB_op": _finite_or_none(db_op),
            "dB_kp": _finite_or_none(db_kp),
        },
        "coeffs": _coeff_block(lin_e),
        "toch": lin_e.toch,
    }
    if problem == "diffraction":
        cs = lin_e.cross_sections()
        block["cross_sections"] = {k: (float(v) if np.isfinite(v) else None) for k, v in cs.items()}
    return block


def _sweep_block(a, eps, miy, sweep_cfg, k0_current):
    if not sweep_cfg or not sweep_cfg.get("enable", True):
        return None
    try:
        k_min = float(sweep_cfg.get("k0_min", 0.25))
        k_max = float(sweep_cfg.get("k0_max", max(2.0 * k0_current, 5.0)))
        npts = int(sweep_cfg.get("points", 160))
    except (TypeError, ValueError):
        raise ComputeError("Свип: некорректные параметры.")
    if not (0 < k_min < k_max and math.isfinite(k_max)):
        raise ComputeError("Свип: требуется 0 < k0_min < k0_max.")
    npts = max(8, min(npts, MAX_SWEEP_POINTS))
    grid = np.linspace(k_min, k_max, npts)
    out = {"k0": [float(x) for x in grid], "q_sca": [], "q_ext": [], "q_abs": [], "q_back": []}
    for k in grid:
        try:
            cs = MieSphere(k0=float(k), a=list(a), eps=list(eps), miy=list(miy),
                           toch=None, problem="diffraction").cross_sections()
        except Exception:
            cs = {}
        for key in ("q_sca", "q_ext", "q_abs", "q_back"):
            v = cs.get(key)
            out[key].append(float(v) if v is not None and np.isfinite(v) else None)
    return out


def compute(payload: dict) -> dict:
    """Главная точка входа: payload (см. app.js) -> словарь для JSON-ответа."""
    try:
        k0 = float(payload.get("k0", 0))
    except (TypeError, ValueError):
        raise ComputeError("k0: некорректное значение.")
    if not (k0 > 0 and math.isfinite(k0)):
        raise ComputeError("k0 должен быть > 0.")

    a, eps, miy = _parse_layers(payload)

    toch = payload.get("toch")
    if toch is not None:
        try:
            toch = int(toch)
        except (TypeError, ValueError):
            raise ComputeError("toch: некорректное значение.")
        if not (1 <= toch <= MAX_TOCH):
            raise ComputeError(f"toch: допустимо 1..{MAX_TOCH}.")

    try:
        phi_deg = float(payload.get("phi_deg", 45.0))
    except (TypeError, ValueError):
        raise ComputeError("phi: некорректное значение.")
    phi_rad = math.radians(phi_deg)

    problems = payload.get("problems") or ["diffraction"]
    problems = [p for p in problems if p in ("diffraction", "antenna")]
    if not problems:
        problems = ["diffraction"]

    results = {}
    for p in problems:
        results[p] = _problem_block(k0, a, eps, miy, toch, phi_rad, p)

    sweep = _sweep_block(a, eps, miy, payload.get("sweep"), k0) \
        if "diffraction" in problems else None

    warnings = []
    sing = abs(math.sin(math.radians(90)) * math.cos(phi_rad))
    if sing > 0.999:
        warnings.append("φ близко к 0°/180°: нормировочный множитель линейной ДН "
                        "сингулярен при θ=90°/270° (особенность оригинальной формулы).")

    return {
        "theta_deg": [float(x) for x in THETA_DEG],
        "results": results,
        "sweep": sweep,
        "meta": {
            "k0": k0,
            "phi_deg": phi_deg,
            "toch_used": results[problems[0]]["toch"],
            "layers_norm": {
                "a": [float(x) for x in a],
                "eps_re": [float(e.real) for e in eps],
                "eps_im": [float(e.imag) for e in eps],
                "mu_re": [float(m.real) for m in miy],
            },
            "warnings": warnings,
        },
    }
