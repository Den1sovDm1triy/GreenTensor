# SPDX-License-Identifier: MIT
"""Сверка ближней зоны по КОМПЛЕКСНЫМ компонентам поля (экспорт hfss_nearvec_*).

Прежний экспорт модуля |E| давал мгновенное значение поля, а не амплитуду; здесь
модуль собирается из Re/Im составляющих: |E| = (|Ex|^2+|Ey|^2+|Ez|^2)^{1/2}.
"""
from __future__ import annotations
import os, sys, math, json
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.abspath(os.path.join("..", "..")))
from green_tensor.nearfield import external_field_sphere as extfield

HERE = os.path.dirname(os.path.abspath(__file__))
HFD = os.path.join(HERE, "reference", "csv")
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)

CLIGHT, F0 = 299792458.0, 10e9
K0 = 2 * math.pi * F0 / CLIGHT
PEC = 1e7
KHAT = (0.0, 0.0, -1.0)

BODIES = {
    "pec":      dict(eps=[PEC],                a_norm=[1.0],                    title="Проводящая сфера"),
    "radome":   dict(eps=[PEC, 2.1 - 4e-4j],   a_norm=[0.9, 1.0],               title="Сфера с укрытием"),
    "luneburg": dict(eps=[1.96, 1.84, 1.64, 1.36, 1.00],
                     a_norm=[0.2, 0.4, 0.6, 0.8, 1.0],                          title="Линза Люнеберга"),
}
CASES = [("pec", "0p5", 0.5), ("pec", "2", 2.0), ("pec", "15", 15.0),
         ("radome", "0p5", 0.5), ("radome", "2", 2.0), ("radome", "15", 15.0),
         ("luneburg", "15", 15.0)]


def db(x):
    return 20.0 * np.log10(np.maximum(np.asarray(x, float), 1e-300))


rows, maps = [], {}
for bk, sk, k0a in CASES:
    tag = f"{bk}_k{sk}_diffraction"
    d = np.loadtxt(os.path.join(HFD, f"hfss_nearvec_{tag}.csv"), delimiter=",", skiprows=1)
    pts_mm, comp = d[:, :3], d[:, 3:]
    E_hf = np.abs(comp[:, 0::2] + 1j * comp[:, 1::2])
    hf = np.linalg.norm(E_hf, axis=1)                    # |E| из Re/Im, В/м

    R = k0a / K0
    b = BODIES[bk]
    Egt = extfield(R, b["eps"], K0, pts_mm * 1e-3, problem="diffraction",
                   a_norm=b["a_norm"], khat=KHAT, total=True, e0=1.0)
    Egt_s = extfield(R, b["eps"], K0, pts_mm * 1e-3, problem="diffraction",
                     a_norm=b["a_norm"], khat=KHAT, total=False, e0=1.0)
    gt_tot = np.linalg.norm(np.abs(Egt), axis=1)
    gt_sca = np.linalg.norm(np.abs(Egt_s), axis=1)

    # вне тела и там, где обе величины значимы
    rr = np.linalg.norm(pts_mm * 1e-3, axis=1)
    ok = (rr > R * 1.02) & np.isfinite(gt_tot) & np.isfinite(gt_sca) & (hf > 0)

    res = {}
    for name, gt in (("полное", gt_tot), ("рассеянное", gt_sca)):
        m = ok & (gt > 0.01 * np.nanmax(gt[ok]))
        res[name] = (float(np.nanmedian(np.abs(db(hf[m]) - db(gt[m])))),
                     float(np.nanmedian(hf[m] / gt[m])),
                     float(db(np.nanmax(hf[ok])) - db(np.nanmax(gt[ok]))))
    best = min(res, key=lambda k: res[k][0])
    gt_best = gt_tot if best == "полное" else gt_sca
    rows.append(dict(tag=tag, body=b["title"], k0a=k0a, kind=best,
                     median_dB=res[best][0], ratio=res[best][1], dpeak_dB=res[best][2],
                     other_kind=[k for k in res if k != best][0],
                     other_median_dB=res[[k for k in res if k != best][0]][0],
                     n_pts=int(ok.sum())))
    maps[tag] = (pts_mm, gt_best, hf, R, b["title"], k0a, best)
    print(f"{tag:26s} {best:11s} медиана={res[best][0]:5.2f} дБ  отношение={res[best][1]:.3f}  "
          f"Δмакс={res[best][2]:+5.2f} дБ  (альтернатива {[k for k in res if k!=best][0]}: "
          f"{res[[k for k in res if k!=best][0]][0]:.2f} дБ)")

with open(os.path.join(HERE, "data", "near_metrics.json"), "w", encoding="utf-8") as fh:
    json.dump(rows, fh, ensure_ascii=False, indent=2)

# --- рисунок: по строке на конфигурацию, слева ТФГ, справа FEM, общая шкала ---
sel = ["pec_k15_diffraction", "radome_k15_diffraction", "luneburg_k15_diffraction"]
fig, axes = plt.subplots(len(sel), 2, figsize=(9.0, 3.6 * len(sel)))
for r, tag in enumerate(sel):
    pts, gt, hf, R, title, k0a, kind = maps[tag]
    ux, uz = np.unique(pts[:, 0]), np.unique(pts[:, 2])
    ix = np.searchsorted(ux, pts[:, 0]); iz = np.searchsorted(uz, pts[:, 2])
    def grid(v):
        g = np.full((uz.size, ux.size), np.nan)
        g[iz, ix] = v
        return np.ma.masked_invalid(g)
    G, H = grid(gt), grid(hf)
    top = float(max(np.nanmax(gt[np.isfinite(gt)]), np.nanmax(hf)))
    nrm = LogNorm(vmin=top * 1e-3, vmax=top)
    cm = plt.get_cmap("inferno").copy(); cm.set_bad("white")
    for c, (M, lab) in enumerate(((G, "метод ТФГ"), (H, "Ansys HFSS (Re/Im)"))):
        ax = axes[r, c]
        im = ax.pcolormesh(ux, uz, M, norm=nrm, cmap=cm, shading="auto")
        ax.add_patch(plt.Circle((0, 0), R * 1e3, fill=False, ls="--", lw=1.0, color="cyan"))
        ax.set_aspect("equal"); ax.set_xlabel("x, мм"); ax.set_ylabel("z, мм")
        ax.set_title(f"{title}, $k_0a={k0a:g}$ — {lab}", fontsize=8.5)
        fig.colorbar(im, ax=ax, label="$|E^s|$, В/м", fraction=0.046)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "near_cmp_vec.png"), dpi=120, bbox_inches="tight")
print("\nрисунок ->", os.path.join(OUT, "near_cmp_vec.png"))
