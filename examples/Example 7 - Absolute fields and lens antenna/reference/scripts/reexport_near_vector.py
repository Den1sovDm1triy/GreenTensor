# SPDX-License-Identifier: MIT
"""Re-export near-field E as VECTOR components (Ex,Ey,Ez, complex) on the x0z grid, from the
already-SOLVED near projects (no re-solve), using the HFSS field calculator directly
(E -> ScalarX/Y/Z -> Real/Imag -> ExportOnGrid). Replaces the scalar Mag_E maps.

Grid: x0z plane (y=0), +/-(a+2*lambda0), NEAR_N points/side (same grid as the Mag_E maps).
Output hfss_nearvec_<tag>.csv: x_mm,y_mm,z_mm, ExRe,ExIm, EyRe,EyIm, EzRe,EzIm  [V/m].

  python reexport_near_vector.py <key> <mode>
  python reexport_near_vector.py all
Diffraction reads EX7_<tag>_near.aedt; antenna reads EX7A_<key>.aedt (Huygens weights re-applied).
"""
import os
import sys

import numpy as np

import run_ex7 as ex7
from ansys.aedt.core import Hfss

try:
    import run_ex7_antenna as ra
except Exception:
    ra = None

COMPS = [("ScalarX", "Real"), ("ScalarX", "Imag"), ("ScalarY", "Real"),
         ("ScalarY", "Imag"), ("ScalarZ", "Real"), ("ScalarZ", "Imag")]
HDR = "x_mm,y_mm,z_mm,ExRe,ExIm,EyRe,EyIm,EzRe,EzIm"


def _export_component(hfss, comp, part, sol, gstart, gstop, gstep, fp):
    fr = hfss.ofieldsreporter
    fr.CalcStack("clear")
    fr.EnterQty("E")
    fr.CalcOp(comp)
    fr.CalcOp(part)
    u = hfss.modeler.model_units
    gsw = [f"{i}{u}" for i in gstart]
    gpw = [f"{i}{u}" for i in gstop]
    gstw = [f"{i}{u}" for i in gstep]
    variation = ["Freq:=", f"{ex7.F_GHZ}GHz"]
    export_options = ["NAME:ExportOption", "IncludePtInOutput:=", True, "RefCSName:=", "Global",
                      "PtInSI:=", True, "FieldInRefCS:=", True]
    fr.ExportOnGrid(fp, gsw, gpw, gstw, sol, variation, export_options,
                    "Cartesian", ["0mm", "0mm", "0mm"], False)


def _col(fp):
    """Return (coords[N,3] mm, value[N]) from a scalar .fld (x,y,z,val in SI metres)."""
    xyz, val = [], []
    for line in open(fp, encoding="utf-8", errors="ignore"):
        p = line.split()
        if len(p) < 4:
            continue
        try:
            v = [float(x) for x in p[:4]]
        except ValueError:
            continue
        xyz.append(v[:3])
        val.append(v[3])
    return np.asarray(xyz) * 1e3, np.asarray(val)


def reexport(key, mode, non_graphical=True):
    v = next(x for x in ex7.CFG["variants"] if x["key"] == key)
    tag = f"{key}_{mode}"
    proj = f"EX7_{tag}_near" if mode == "diffraction" else f"EX7A_{key}"
    ppath = os.path.join(ex7.OUT, proj + ".aedt")
    if not os.path.exists(ppath):
        print(f"[skip] {tag}: project not found ({proj})")
        return False
    hfss = Hfss(project=ppath, design=proj, version=ex7.AEDT_VERSION,
                non_graphical=non_graphical, new_desktop=True)
    try:
        sol = "Setup1 : LastAdaptive"
        if mode == "antenna" and ra is not None:
            try:
                scale, c = ra.load_calib()
                ra._set(hfss, "H_Ie", ["H_Im"], scale, scale * c)
            except Exception as e:
                print(f"   [warn] source weights: {e}")
        h = v["near_half_mm"]
        step = 2 * h / (ex7.NEAR_N - 1)
        cols, coords = [], None
        for comp, part in COMPS:
            fp = os.path.join(ex7.OUT, f"_c_{tag}_{comp}_{part}.fld")
            _export_component(hfss, comp, part, sol, [-h, 0.0, -h], [h, 0.0, h], [step, step, step], fp)
            xyz, val = _col(fp)
            if coords is None:
                coords = xyz
            cols.append(val)
            os.remove(fp)
        out = np.column_stack([coords] + cols)
        csvp = os.path.join(ex7.OUT, f"hfss_nearvec_{tag}.csv")
        with open(csvp, "w", encoding="utf-8", newline="") as f:
            f.write(HDR + "\n")
            for row in out:
                f.write(",".join(f"{x:.6g}" for x in row) + "\n")
        print(f"[ok] {tag}: hfss_nearvec_{tag}.csv ({out.shape[0]} pts, {out.shape[1]} cols)")
        return True
    finally:
        try:
            hfss.release_desktop(close_projects=False, close_desktop=True)
        except Exception:
            pass


if __name__ == "__main__":
    if sys.argv[1] == "all":
        keys = [v["key"] for v in ex7.CFG["variants"]]
        jobs = [(k, "diffraction") for k in keys]
        jobs += [(k, "antenna") for k in ("lens_k0p5", "lens_k2", "lens_k15", "luneburg_k15")]
        for k, m in jobs:
            try:
                reexport(k, m)
            except Exception as e:
                print(f"[FAIL] {k}_{m}: {type(e).__name__}: {str(e)[:150]}")
    else:
        reexport(sys.argv[1], sys.argv[2])
    print("Done.")
