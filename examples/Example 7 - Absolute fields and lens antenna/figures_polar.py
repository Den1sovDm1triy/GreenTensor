# SPDX-License-Identifier: MIT
"""Круговые (полярные) диаграммы: ТФГ против полноволнового эталона Ansys HFSS.

Для каждой конфигурации строится набор из четырёх диаграмм 2x2:
  верхний ряд  -- линейная поляризация: вертикальная (E-плоскость, phi=0)
                  и горизонтальная (H-плоскость, phi=90 град);
  нижний  ряд  -- круговая поляризация: основная (ОП) и кросс-поляризация (КП).
На каждой диаграмме -- одна пара кривых: расчёт ТФГ и сопоставимый эталон HFSS.
Ориентация всех диаграмм одинакова: theta=0 (обратное рассеяние) справа,
отсчёт против часовой стрелки, theta=180 град (рассеяние вперёд) слева.

HFSS отдаёт rE [В] на сфере r = R_OBS -> |E| = rE/R_OBS.
"""
from __future__ import annotations
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
GTD = os.path.join(HERE, "data")
HFD = os.path.join(HERE, "reference", "csv")
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
R_OBS = 10.0
DYN_DB = 45.0          # показываемый динамический диапазон полярной диаграммы

C_GT, C_HF = "#b3202c", "#1f5fa8"

CASES = [
    # (тег, подпись тела, k0a, режим)
    ("pec_k0p5_diffraction",      "Проводящая сфера",          0.5,  "diffraction"),
    ("pec_k2_diffraction",        "Проводящая сфера",          2.0,  "diffraction"),
    ("pec_k15_diffraction",       "Проводящая сфера",          15.0, "diffraction"),
    ("radome_k0p5_diffraction",   "Сфера с укрытием",          0.5,  "diffraction"),
    ("radome_k2_diffraction",     "Сфера с укрытием",          2.0,  "diffraction"),
    ("radome_k15_diffraction",    "Сфера с укрытием",          15.0, "diffraction"),
    ("luneburg_k15_diffraction",  "Линза Люнеберга",           15.0, "diffraction"),
]


def db(x):
    return 20.0 * np.log10(np.maximum(np.asarray(x, float), 1e-300))


def load(tag):
    g = np.loadtxt(os.path.join(GTD, f"far_{tag}.csv"), delimiter=",", skiprows=1)
    h = np.loadtxt(os.path.join(HFD, f"hfss_far_{tag}.csv"), delimiter=",", skiprows=1)
    c = np.loadtxt(os.path.join(HFD, f"hfss_circ_{tag}.csv"), delimiter=",", skiprows=1)
    return dict(
        gth=g[:, 0], g_v=g[:, 1], g_h=g[:, 2], g_x=g[:, 3], g_op=g[:, 4], g_kp=g[:, 5],
        hth=h[:, 0], h_v=h[:, 1] / R_OBS, h_h=h[:, 2] / R_OBS, h_x=h[:, 3] / R_OBS,
        cth=c[:, 0], h_lhcp=c[:, 1] / R_OBS, h_rhcp=c[:, 2] / R_OBS)


def metrics(gth, gv, hth, hv, floor=0.01):
    """Медиана |Δ| в дБ по значимым лепесткам и расхождение максимумов."""
    hi = np.interp(gth % 360, hth % 360, hv, period=360)
    m = gv > floor * np.nanmax(gv)
    return (float(np.nanmedian(np.abs(db(hi[m]) - db(gv[m])))),
            float(db(np.nanmax(hi)) - float(db(np.nanmax(gv)))))


def panel(ax, th_g, v_g, th_h, v_h, title, top_db):
    lo = top_db - DYN_DB
    ax.plot(np.deg2rad(th_g), np.clip(db(v_g), lo, None), lw=2.1, color=C_GT, label="ТФГ")
    ax.plot(np.deg2rad(th_h), np.clip(db(v_h), lo, None), lw=1.7, ls=(0, (5, 2)),
            color=C_HF, label="Ansys HFSS")
    ax.set_theta_zero_location("E"); ax.set_theta_direction(1)
    ax.set_rlim(lo, top_db)
    ax.set_rlabel_position(78.0)
    ax.set_rticks(np.arange(lo + 15.0, top_db - 0.1, 15.0))
    ax.set_thetagrids(np.arange(0, 360, 30),
                      [f"{t}°" if t in (0, 180) else "" for t in range(0, 360, 30)],
                      fontsize=15)
    ax.tick_params(axis="y", labelsize=11.5, pad=-2)
    ax.grid(alpha=0.4, lw=0.5)
    ax.set_title(title, fontsize=16, pad=9)


report = {}
for tag, body, k0a, prob in CASES:
    d = load(tag)
    # рукость круговой поляризации: фиксированное соответствие проверяется
    # сопоставлением обоих вариантов, выигравший записывается в отчёт
    a1 = 0.5 * (metrics(d["gth"], d["g_op"], d["cth"], d["h_rhcp"])[0]
                + metrics(d["gth"], d["g_kp"], d["cth"], d["h_lhcp"])[0])
    a2 = 0.5 * (metrics(d["gth"], d["g_op"], d["cth"], d["h_lhcp"])[0]
                + metrics(d["gth"], d["g_kp"], d["cth"], d["h_rhcp"])[0])
    op_is_rhcp = a1 <= a2
    h_op = d["h_rhcp"] if op_is_rhcp else d["h_lhcp"]
    h_kp = d["h_lhcp"] if op_is_rhcp else d["h_rhcp"]

    pairs = [
        ("вертикальная, $\\varphi=0$",            d["g_v"],  d["h_v"]),
        ("горизонтальная, $\\varphi=90^\\circ$",   d["g_h"],  d["h_h"]),
        ("круговая, основная (ОП)",             d["g_op"], h_op),
        ("круговая, кросс (КП)",                d["g_kp"], h_kp),
    ]

    # единый верхний предел шкалы на всю четвёрку -> диаграммы сопоставимы
    top = np.ceil(max(db(np.nanmax(v)) for _, v, _ in pairs) / 5.0) * 5.0

    fig, axes = plt.subplots(2, 2, figsize=(8.6, 7.7),
                             subplot_kw=dict(projection="polar"))
    met = {}
    for ax, (name, gv, hv) in zip(axes.ravel(), pairs):
        th_h = d["cth"] if name.startswith("круговая") else d["hth"]
        # в антенном режиме нормировка источника расходится -> сверяется ФОРМА,
        # абсолютное смещение уровня выносится отдельной величиной
        md, dpk = metrics(d["gth"], gv, th_h, hv)
        key = ("v", "h", "op", "kp")[[p[0] for p in pairs].index(name)]
        met[key] = dict(median_dB=md, dpeak_dB=dpk,
                        gt_max_Vm=float(np.nanmax(gv)), hf_max_Vm=float(np.nanmax(hv)))
        panel(ax, d["gth"], gv, th_h, hv,
              f"{name}\n$\\Delta$пик {dpk:+.2f} дБ, медиана {md:.2f} дБ", top)
    hs, ls = axes[0, 0].get_legend_handles_labels()
    fig.legend(hs, ls, loc="lower center", ncol=2, fontsize=15, frameon=False,
               bbox_to_anchor=(0.5, 0.005))
    # заголовок набора не рисуется: тело, размер, режим и шкала указаны в подписи
    fig.subplots_adjust(left=0.05, right=0.97, top=0.96, bottom=0.085,
                        wspace=0.20, hspace=0.44)
    fig.savefig(os.path.join(OUT, f"pol_{tag}.png"), dpi=125, bbox_inches="tight")
    plt.close(fig)

    met["xpol_gt_Vm"] = float(np.nanmax(d["g_x"]))
    met["xpol_hf_Vm"] = float(np.nanmax(d["h_x"]))
    met["op_is_rhcp"] = bool(op_is_rhcp)
    met["body"], met["k0a"], met["problem"] = body, k0a, prob
    report[tag] = met
    print(f"{tag:26s} op↔{'RHCP' if op_is_rhcp else 'LHCP':4s} "
          + "  ".join(f"{k}:{met[k]['median_dB']:5.2f}/{met[k]['dpeak_dB']:+5.2f}"
                      for k in ("v", "h", "op", "kp")))

with open(os.path.join(GTD, "polar_metrics.json"), "w", encoding="utf-8") as fh:
    json.dump(report, fh, ensure_ascii=False, indent=2)
print("\nfigures ->", OUT)
