# SPDX-License-Identifier: MIT
"""Антенный режим линзы Люнеберга: ДН по теореме взаимности против эталона.

Срезы E и H, шкала АБСОЛЮТНАЯ (В/м); расчёт --
``green_tensor.antenna_pattern_sphere``.
"""
from __future__ import annotations
import os, sys, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))
from green_tensor.nearfield import antenna_pattern_sphere

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
HFD = os.path.join(HERE, "reference", "csv")
DYN, C_GT, C_HF = 45.0, "#b3202c", "#1f5fa8"

a = 71.5702e-3
K0 = 2 * math.pi * 10e9 / 299792458.0
LAM = 2 * math.pi / K0
EPS = [1.96, 1.84, 1.64, 1.36, 1.00]
ANORM = [0.2, 0.4, 0.6, 0.8, 1.0]
GAP = 0.074948115e-3

h = np.loadtxt(os.path.join(HFD, "hfss_far_luneburg_k15_antenna.csv"),
               delimiter=",", skiprows=1)
hth, hE, hH = h[:, 0], h[:, 1] / 10.0, h[:, 2] / 10.0

th = np.deg2rad(np.arange(0.25, 360.0, 0.25))
thd = np.degrees(th)
rec = {p: antenna_pattern_sphere(a, EPS, K0, th, current=1.0, length=LAM / 50.0,
                                 gap=GAP, r=10.0, plane=p, a_norm=ANORM)["E"]
       for p in ("E", "H")}


def db(x):
    return 20.0 * np.log10(np.maximum(np.asarray(x, float), 1e-300))


def metrics(td, v, rth, rv, norm=False):
    hi = np.interp(td % 360, rth % 360, rv, period=360)
    m = v > 0.01 * np.nanmax(v)
    x, y = (v / np.nanmax(v), hi / np.nanmax(hi)) if norm else (v, hi)
    return (float(np.nanmedian(np.abs(db(y[m]) - db(x[m])))),
            float(db(np.nanmax(hi)) - db(np.nanmax(v))))


def panel(ax, td, v, rth, rv, title, top, norm):
    a_, b_ = (v / np.nanmax(v), rv / np.nanmax(rv)) if norm else (v, rv)
    lo = top - DYN
    ax.plot(np.deg2rad(td), np.clip(db(a_), lo, None), lw=2.1, color=C_GT, label="ТФГ")
    ax.plot(np.deg2rad(rth), np.clip(db(b_), lo, None), lw=1.7, ls=(0, (5, 2)),
            color=C_HF, label="Ansys HFSS")
    ax.set_theta_zero_location("E"); ax.set_theta_direction(1)
    ax.set_rlim(lo, top); ax.set_rlabel_position(78.0)
    ax.set_rticks(np.arange(lo + 15.0, top - 0.1, 15.0))
    ax.set_thetagrids(np.arange(0, 360, 30),
                      [f"{t}°" if t in (0, 180) else "" for t in range(0, 360, 30)],
                      fontsize=15)
    ax.tick_params(axis="y", labelsize=11.5, pad=-2)
    ax.grid(alpha=0.4, lw=0.5); ax.set_title(title, fontsize=16, pad=9)


fig, ax = plt.subplots(1, 2, figsize=(8.6, 4.3), subplot_kw=dict(projection="polar"))
top = np.ceil(db(max(np.nanmax(rec["E"]), np.nanmax(hE))) / 5.0) * 5.0
for k, (pl, rv, nm) in enumerate((("E", hE, "вертикальная, $\\varphi=0$"),
                                  ("H", hH, "горизонтальная, $\\varphi=90^\\circ$"))):
    med, dpk = metrics(thd, rec[pl], hth, rv)
    panel(ax[k], thd, rec[pl], hth, rv,
          f"{nm}\n$\\Delta$пик {dpk:+.2f} дБ, медиана {med:.2f} дБ", top, False)
hs, ls = ax[0].get_legend_handles_labels()
fig.legend(hs, ls, loc="lower center", ncol=2, fontsize=15, frameon=False,
           bbox_to_anchor=(0.5, 0.0))
fig.subplots_adjust(left=0.05, right=0.97, top=0.90, bottom=0.14, wspace=0.20)
fig.savefig(os.path.join(OUT, "pol_luneburg_k15_antenna.png"), dpi=125,
            bbox_inches="tight")

for pl, rv in (("E", hE), ("H", hH)):
    med, dpk = metrics(thd, rec[pl], hth, rv)
    print(f"  взаимность, {pl}-плоскость: Дпик={dpk:+.2f} дБ, медиана={med:.2f} дБ")
print("рисунок ->", os.path.join(OUT, "pol_luneburg_k15_antenna.png"))
