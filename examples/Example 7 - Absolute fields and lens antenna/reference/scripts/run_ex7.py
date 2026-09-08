# SPDX-License-Identifier: MIT
"""Example 7 (see ../MODELS.md) — absolute fields, FEM (Driven Modal), PyAEDT 1.4 / AEDT 2024.2.

3 bodies (pec / radome=PEC+PTFE / lens eps=4) x 3 sizes (k0a=0.5,2,15) x 2 modes
(diffraction = plane wave E0=1 V/m; antenna = elementary current dipole I=1 A on the
surface pole). Everything ABSOLUTE (V, V/m) — no peak normalisation anywhere.

Per (body,size,mode) tag = <body>_k{0p5|2|15}_{diffraction|antenna} it writes into EX7_KIT/:
  hfss_far_<tag>.csv     Theta, rETheta(phi0), rEPhi(phi90), rETheta(phi90=cross)   [V]
  hfss_circ_<tag>.csv    Theta, rELHCP, rERHCP  (phi0 cut)                           [V]
  hfss_near_<tag>.csv    x,y,z, Mag_E on the x0z plane grid                          [V/m]
  hfss_farmap_<tag>.csv  Theta, Phi, rETotal   (theta 0..180/3, phi 0..360/5)       [V]
  hfss_profile_<tag>.txt peak memory / wall time / tetrahedra
A per-tag status dict is collected into EX7_KIT/ex7_run_summary.json.

rE is scaled to SI (r in metres): rE_SI = rE_HFSS(mm)*1e-3, so |E|(r=10m) = rE_SI/10.
Circular is SYNTHESISED from two orthogonal linear solves (x- and y-pol) for diffraction,
and from crossed x/y dipoles with 90 deg phase for the antenna (native rELHCP/rERHCP).

Run one:   python run_ex7.py pec_k0p5 diffraction
Run body:  python run_ex7.py pec            # all sizes+modes for a body
Run all:   python run_ex7.py
Flags:     --graphical  --no-circ  --no-maps  --modes diffraction,antenna
"""
import argparse
import glob
import json
import math
import os
import re
import shutil
import sys
import time

import numpy as np

try:  # console/redirect may be cp1251 on RU Windows -> Greek eps in titles crashes print()
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "EX7_KIT")
CFG = json.load(open(os.path.join(HERE, "ex7_variants.json"), encoding="utf-8"))
AEDT_VERSION = os.environ.get("AEDT_VERSION", "2024.2")
CORES = int(os.environ.get("HFSS_CORES", "16"))   # 7950X3D: 16 physical cores; license unlimited
MAXPASSES = int(os.environ.get("HFSS_MAXPASSES", "12"))
MAXDELTAS = float(os.environ.get("HFSS_MAXDELTAS", "0.02"))
NEAR_N = int(os.environ.get("EX7_NEAR_N", "101"))          # near-field grid points/side
RE_SCALE = 1e-3                                            # mm -> m (r factor to SI)
# Antenna elementary-source calibration: (ideal Hertzian I*l moment)/(HFSS current-sheet moment),
# measured in free space (guide sec.6 "domnozhit' na I_treb/I_fakt"). 1.0 = diffraction / uncalibrated.
ANT_CAL = float(os.environ.get("EX7_ANT_CAL", "1.0"))
AMP = 1.0                                                  # per-run amplitude factor (set for antenna)
OBS_R = float(CFG.get("obs_r_m", 10.0))
F_HZ = CFG["frequency_Hz"]
F_GHZ = F_HZ / 1e9
LAM0 = CFG["lambda0_mm"]
DIP_L = CFG["antenna"]["dipole_len_mm"]                    # lambda/50

from ansys.aedt.core import Hfss


# --------------------------------------------------------------------------- geometry
def _is_vacuum(mat):
    return (mat.get("kind") == "dielectric" and float(mat.get("eps_re", 1.0)) == 1.0
            and float(mat.get("tand_e", 0.0)) == 0.0)


def _ensure_material(hfss, mat):
    if mat["kind"] == "pec":
        return "pec"
    name = f"diel_e{mat['eps_re']:g}_t{mat['tand_e']:g}"
    if hfss.materials.exists_material(name):
        return name
    m = hfss.materials.add_material(name)
    m.permittivity = mat["eps_re"]
    m.dielectric_loss_tangent = mat["tand_e"]
    m.update()
    return name


def build_body(hfss, v):
    built = []
    for lyr in v["layers_hfss"]:
        mat = lyr["material"]
        if _is_vacuum(mat):
            continue
        matname = _ensure_material(hfss, mat)
        obj = hfss.modeler.create_sphere([0, 0, 0], lyr["r_outer_mm"], name=f"L{lyr['index']}",
                                         material=matname)
        if lyr["r_inner_mm"] and lyr["r_inner_mm"] > 0:
            cut = hfss.modeler.create_sphere([0, 0, 0], lyr["r_inner_mm"], name=f"cut{lyr['index']}")
            hfss.modeler.subtract(obj, cut, keep_originals=False)
            obj.material_name = matname
        built.append(f"L{lyr['index']}")
    return built


# --------------------------------------------------------------------------- dipole
def build_dipole(hfss, axis, a, name):
    """Elementary current element (I=1 A, length lambda/50) at the +z pole surface.
    axis 'z' = radial (vertical); 'x','y' = tangential. Small standoff g avoids PEC short."""
    l, w = DIP_L, DIP_L / 8.0
    g = DIP_L * 0.25
    z0 = a + g
    if axis == "z":
        sh = hfss.modeler.create_rectangle("XZ", [-w / 2, 0.0, z0], [l, w], name=f"dip_{name}")
        il = [[0.0, 0.0, z0], [0.0, 0.0, z0 + l]]
    elif axis == "x":
        sh = hfss.modeler.create_rectangle("XY", [-l / 2, -w / 2, z0], [l, w], name=f"dip_{name}")
        il = [[-l / 2, 0.0, z0], [l / 2, 0.0, z0]]
    else:  # y
        sh = hfss.modeler.create_rectangle("XY", [-w / 2, -l / 2, z0], [w, l], name=f"dip_{name}")
        il = [[0.0, -l / 2, z0], [0.0, l / 2, z0]]
    src = hfss.assign_current_source_to_sheet(sh.name, start_direction=il, name=f"I_{name}")
    return src


# --------------------------------------------------------------------------- far field
def _far_complex(hfss, sol, phi):
    """Complex (rETheta, rEPhi) [SI V] and theta[deg] for one phi cut of setup 'ff'."""
    data = hfss.post.get_solution_data(
        expressions=["rETheta", "rEPhi"], setup_sweep_name=sol, primary_sweep_variable="Theta",
        variations={"Theta": ["All"], "Phi": [f"{phi}deg"]},
        report_category="Far Fields", context="ff")
    if not data:
        raise RuntimeError(f"no far-field data at phi={phi}")
    theta = np.asarray(data.primary_sweep_values, dtype=float)

    def cplx(expr):
        r = np.asarray(data.get_expression_data(expr, formula="real")[1], float).ravel()
        i = np.asarray(data.get_expression_data(expr, formula="imag")[1], float).ravel()
        return (r + 1j * i) * RE_SCALE * AMP
    return theta, cplx("rETheta"), cplx("rEPhi")


def _write(path, header, cols):
    arrs = [np.asarray(c, float) for c in cols]
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(header + "\n")
        for row in zip(*arrs):
            f.write(",".join(f"{x:.9g}" for x in row) + "\n")


def export_far(hfss, sol, tag):
    """hfss_far_<tag>.csv : rETheta(phi0), rEPhi(phi90), rETheta(phi90 cross)."""
    th0, et0, _ = _far_complex(hfss, sol, 0)
    th90, et90, ep90 = _far_complex(hfss, sol, 90)
    _write(os.path.join(OUT, f"hfss_far_{tag}.csv"),
           "Theta,rETheta_phi0_V,rEPhi_phi90_V,rETheta_phi90_cross_V",
           [th0, np.abs(et0), np.abs(ep90), np.abs(et90)])
    return th0, et0  # phi0 rETheta for the reference check


def _circ_write(th, etx, epx, ety, epy, tag):
    """Circular-incidence synthesis on phi=0. Incident LHCP=(x - j y)/sqrt2 (guide E=(1,-j,0));
    scattered circular comps rE_R=(Ftheta - j Fphi)/sqrt2, rE_L=(Ftheta + j Fphi)/sqrt2."""
    s2 = math.sqrt(2.0)
    fth = (etx - 1j * ety) / s2
    fph = (epx - 1j * epy) / s2
    reR = np.abs((fth - 1j * fph) / s2)
    reL = np.abs((fth + 1j * fph) / s2)
    _write(os.path.join(OUT, f"hfss_circ_{tag}.csv"), "Theta,rELHCP_V,rERHCP_V", [th, reL, reR])


def export_circ_native(hfss, sol, tag):
    """Antenna circular: HFSS rELHCP/rERHCP directly (crossed dipoles, 90 deg phase)."""
    data = hfss.post.get_solution_data(
        expressions=["rELHCP", "rERHCP"], setup_sweep_name=sol, primary_sweep_variable="Theta",
        variations={"Theta": ["All"], "Phi": ["0deg"]}, report_category="Far Fields", context="ff")
    if not data:
        raise RuntimeError("no circular far-field data")
    th = np.asarray(data.primary_sweep_values, float)
    rl = np.abs(np.asarray(data.get_expression_data("rELHCP", formula="real")[1], float).ravel()
                + 1j * np.asarray(data.get_expression_data("rELHCP", formula="imag")[1], float).ravel()) * RE_SCALE * AMP
    rr = np.abs(np.asarray(data.get_expression_data("rERHCP", formula="real")[1], float).ravel()
                + 1j * np.asarray(data.get_expression_data("rERHCP", formula="imag")[1], float).ravel()) * RE_SCALE * AMP
    _write(os.path.join(OUT, f"hfss_circ_{tag}.csv"), "Theta,rELHCP_V,rERHCP_V", [th, rl, rr])


def export_farmap(hfss, sol, tag):
    """hfss_farmap_<tag>.csv : rETotal(theta,phi) on the map sphere 'ffmap' (theta 0..180/3,
    phi 0..360/5). Iterate phi explicitly for a clean (Theta,Phi,rETotal) table."""
    path = os.path.join(OUT, f"hfss_farmap_{tag}.csv")
    wrote = 0
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("Theta,Phi,rETotal_V\n")
        for p in range(0, 361, 5):
            d2 = hfss.post.get_solution_data(
                expressions=["rETotal"], setup_sweep_name=sol, primary_sweep_variable="Theta",
                variations={"Theta": ["All"], "Phi": [f"{p}deg"]},
                report_category="Far Fields", context="ffmap")
            if not d2:
                continue
            tt = np.asarray(d2.primary_sweep_values, float)
            yy = np.asarray(d2.get_expression_data("rETotal", formula="real")[1], float).ravel() * RE_SCALE * AMP
            for t, y in zip(tt, yy):
                f.write(f"{t:.6g},{p:.6g},{y:.9g}\n")
                wrote += 1
    if not wrote:
        raise RuntimeError("farmap produced no rows")


def export_near_vector(hfss, sol, tag, v):
    """hfss_nearvec_<tag>.csv : COMPLEX Ex,Ey,Ez on the x0z grid (V/m). Uses the field calculator
    (E -> ScalarX/Y/Z -> Real/Imag) as named expressions, exported on the grid IN-SESSION (the
    volumetric field is not saved to disk, so this must run right after the solve)."""
    h = v["near_half_mm"]
    step = 2 * h / (NEAR_N - 1)
    fr = hfss.ofieldsreporter
    specs = [("nvExRe", "ScalarX", "Real"), ("nvExIm", "ScalarX", "Imag"),
             ("nvEyRe", "ScalarY", "Real"), ("nvEyIm", "ScalarY", "Imag"),
             ("nvEzRe", "ScalarZ", "Real"), ("nvEzIm", "ScalarZ", "Imag")]
    coords, cols = None, []
    for name, comp, part in specs:
        fr.CalcStack("clear")
        fr.EnterQty("E")
        fr.CalcOp(comp)
        fr.CalcOp(part)
        try:
            fr.AddNamedExpression(name, "Fields")
        except Exception:
            pass
        fld = os.path.join(OUT, f"_nv_{tag}_{name}.fld")
        hfss.post.export_field_file_on_grid(
            name, solution=sol, variations={}, file_name=fld, grid_type="Cartesian",
            grid_start=[-h, 0.0, -h], grid_stop=[h, 0.0, h], grid_step=[step, step, step],
            is_vector=False, intrinsics={"Freq": f"{F_GHZ}GHz"},
            export_with_sample_points=True, export_in_si_system=True)
        xyz, val = [], []
        for line in open(fld, encoding="utf-8", errors="ignore"):
            p = line.split()
            if len(p) >= 4:
                try:
                    q = [float(x) for x in p[:4]]
                except ValueError:
                    continue
                xyz.append(q[:3])
                val.append(q[3])
        if coords is None:
            coords = xyz
        cols.append(val)
        try:
            os.remove(fld)
        except OSError:
            pass
    csv = os.path.join(OUT, f"hfss_nearvec_{tag}.csv")
    with open(csv, "w", encoding="utf-8", newline="") as f:
        f.write("x_mm,y_mm,z_mm,ExRe,ExIm,EyRe,EyIm,EzRe,EzIm\n")
        for i in range(len(coords)):
            x, y, z = coords[i]
            comps = ",".join(f"{cols[k][i]:.6g}" for k in range(6))
            f.write(f"{x*1e3:.6g},{y*1e3:.6g},{z*1e3:.6g},{comps}\n")
    return len(coords)


def export_near(hfss, sol, tag, v):
    """hfss_near_<tag>.csv : Mag_E on the x0z plane, +/-(a+2lambda), NEAR_N points/side.
    Also writes hfss_nearvec_<tag>.csv with complex Ex,Ey,Ez on the same grid."""
    h = v["near_half_mm"]
    step = 2 * h / (NEAR_N - 1)
    f_fld = os.path.join(OUT, f"_near_{tag}.fld")
    hfss.post.export_field_file_on_grid(
        "Mag_E", solution=sol, variations={},
        file_name=f_fld, grid_type="Cartesian",
        grid_start=[-h, 0.0, -h], grid_stop=[h, 0.0, h], grid_step=[step, step, step],
        is_vector=False, intrinsics={"Freq": f"{F_GHZ}GHz", "Phase": "0deg"},
        export_with_sample_points=True, export_in_si_system=True)
    # convert the .fld to a clean CSV (x,y,z,Mag_E)
    csv = os.path.join(OUT, f"hfss_near_{tag}.csv")
    rows = []
    for line in open(f_fld, encoding="utf-8", errors="ignore"):
        parts = line.split()
        if len(parts) >= 4:
            try:
                vals = [float(x) for x in parts[:4]]
                rows.append(vals)
            except ValueError:
                continue
    with open(csv, "w", encoding="utf-8", newline="") as f:
        f.write("x_mm,y_mm,z_mm,Mag_E_Vm\n")
        for r in rows:
            # .fld sample points are SI metres -> convert x,y,z to mm; Mag_E already V/m (x AMP for antenna)
            f.write(f"{r[0]*1e3:.6g},{r[1]*1e3:.6g},{r[2]*1e3:.6g},{r[3]*AMP:.6g}\n")
    try:
        os.remove(f_fld)
    except OSError:
        pass
    try:
        nv = export_near_vector(hfss, sol, tag, v)
        print(f"    [nearvec] {tag}: Ex,Ey,Ez exported ({nv} pts)")
    except Exception as exc:
        print(f"    [nearvec] {tag}: vector export failed: {exc}")
    return len(rows)


# --------------------------------------------------------------------------- profile
def export_profile(hfss, setup, tag, extra):
    info = dict(extra)
    base = os.path.join(OUT, f"hfss_profile_{tag}")
    peak_mb, tets = None, None
    try:
        pf = hfss.export_profile("Setup1", output_file=base + ".prof")
        txt = open(pf, encoding="utf-8", errors="ignore").read()
        mm = 0.0
        for mo in re.finditer(r"[Mm]ax(?:imum)? memory/process\s*[:=]?\s*([\d.]+)\s*([GMK])", txt):
            u = mo.group(2)
            mm = max(mm, float(mo.group(1)) * (1024 if u == "G" else (1 if u == "M" else 1 / 1024)))
        peak_mb = round(mm, 1) if mm else None
        mo = re.search(r"Max(?:imum)? solved tets?\s*[:=]?\s*([\d,]+)", txt)
        if mo:
            tets = int(mo.group(1).replace(",", ""))
    except Exception as e:
        info["profile_err"] = str(e)[:150]
    try:
        ms = hfss.export_mesh_stats("Setup1", output_file=base + ".mstat")
        for line in open(ms, encoding="utf-8", errors="ignore"):
            mo = re.search(r"Total number of mesh elements\s*[:=]?\s*([\d,]+)", line)
            if mo:
                tets = max(tets or 0, int(mo.group(1).replace(",", "")))
    except Exception:
        pass
    info["peak_memory_MB"] = peak_mb
    info["tetrahedra"] = tets
    with open(base + ".txt", "w", encoding="utf-8") as f:
        f.write(f"# HFSS FEM profile — {tag}\n")
        for k, val in info.items():
            f.write(f"{k}: {val}\n")
    return info


# --------------------------------------------------------------------------- setup
def _clean(proj):
    for p in glob.glob(os.path.join(OUT, proj + ".aedt*")) + glob.glob(os.path.join(OUT, proj + ".aedb")):
        try:
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
        except Exception:
            pass
    r = os.path.join(OUT, proj + ".aedtresults")
    if os.path.isdir(r):
        shutil.rmtree(r, ignore_errors=True)


PAD_LAM = float(os.environ.get("EX7_PAD_LAM", "2.0"))         # near-field air pad beyond surface [lambda0]
PAD_SMALL_LAM = float(os.environ.get("EX7_PAD_SMALL", "0.25"))  # far-field open-region pad [lambda0] (~15mm gap)
NEAR_PASSES = int(os.environ.get("EX7_NEAR_PASSES", "6"))     # adaptive passes for the (big) near solve


def _open_region(hfss, pad_lam):
    """Open Region (scattered-field formulation). Auto-pad = quarter-wave at OpFreq, so a lower
    OpFreq enlarges the domain: pad = (pad_lam+0.25)*lambda0 -> OpFreq = c/(4*pad)."""
    pad_m = (pad_lam + 0.25) * LAM0 / 1000.0
    f_op = 299792458.0 / (4.0 * pad_m)
    hfss.create_open_region(frequency=f"{max(f_op, 1e6) / 1e9:.4f}GHz")


def _setup(hfss, passes):
    s = hfss.create_setup("Setup1")
    s.props["Frequency"] = f"{F_GHZ}GHz"
    s.props["MaximumPasses"] = passes
    s.props["MaxDeltaS"] = MAXDELTAS
    s.props["MinimumConvergedPasses"] = 2
    try:
        s.props["BasisOrder"] = 1              # First Order (guide sec.4.4)
    except Exception:
        pass
    s.update()
    return s


def _excite(hfss, v, mode, pol="x"):
    if mode == "diffraction":
        p = [1, 0, 0] if pol == "x" else [0, 1, 0]
        hfss.plane_wave(vector_format="Cartesian", origin=[0, 0, 0], polarization=p,
                        propagation_vector=[0, 0, -1], name=f"pw_{pol}")
    else:
        build_dipole(hfss, "z", v["a_mm"], "z")   # radial (vertical) elementary dipole, I=1 A


def _new_hfss(proj, non_graphical):
    _clean(proj)
    return Hfss(project=os.path.join(OUT, proj + ".aedt"), design=proj, solution_type="Modal",
                version=AEDT_VERSION, non_graphical=non_graphical, new_desktop=True)


# --------------------------------------------------- phase A: far field / circular / far map
def phase_ff(v, mode, tag, non_graphical, do_circ, do_maps, res):
    global AMP
    AMP = ANT_CAL if mode == "antenna" else 1.0
    hfss = _new_hfss(f"EX7_{tag}", non_graphical)
    try:
        hfss.modeler.model_units = "mm"
        build_body(hfss, v)
        a = v["a_mm"]
        _excite(hfss, v, mode, "x")
        _open_region(hfss, PAD_SMALL_LAM)                       # small, accurate scattered FF
        hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=360, theta_step=1,
                                    phi_start=0, phi_stop=90, phi_step=90, name="ff")
        hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=180, theta_step=3,
                                    phi_start=0, phi_stop=360, phi_step=5, name="ffmap")
        _setup(hfss, MAXPASSES)
        hfss.save_project()
        print(f"[i] {tag}: FF solve (cores={CORES}, passes={MAXPASSES}) ...")
        hfss.analyze_setup("Setup1", cores=CORES)
        sol = "Setup1 : LastAdaptive"

        th0, et0 = export_far(hfss, sol, tag)
        res["outputs"].append(f"hfss_far_{tag}.csv")
        ref = v["reference_Vm"][mode]
        pkE = float(np.max(np.abs(et0))) / OBS_R
        res["checks"]["peak_Eplane_Vm"] = pkE
        res["checks"]["ref_Eplane_Vm"] = ref["Eplane"]
        print(f"[check] {tag}: peak |E_theta|(phi0)/10m = {pkE:.4g} V/m  (guide ref {(ref['Eplane'] if ref['Eplane'] is not None else float('nan')):.4g})")

        if do_maps:
            try:
                export_farmap(hfss, sol, tag)
                res["outputs"].append(f"hfss_farmap_{tag}.csv")
            except Exception as e:
                res["checks"]["farmap_err"] = str(e)[:160]
                print(f"[!] {tag}: far map failed: {e}")

        if do_circ:
            try:
                if mode == "diffraction":
                    thc, etx0, epx0 = _far_complex(hfss, sol, 0)
                    for b in hfss.boundaries:
                        if b.name == "pw_x":
                            b.delete()
                    _excite(hfss, v, mode, "y")
                    hfss.analyze_setup("Setup1", cores=CORES)
                    _, ety0, epy0 = _far_complex(hfss, "Setup1 : LastAdaptive", 0)
                    _circ_write(thc, etx0, epx0, ety0, epy0, tag)
                else:
                    build_dipole(hfss, "x", a, "x")
                    build_dipole(hfss, "y", a, "y")
                    hfss.analyze_setup("Setup1", cores=CORES)
                    hfss.edit_sources({"I_x": ("1A", "0deg"), "I_y": ("1A", "90deg"),
                                       "I_z": ("0A", "0deg")})
                    export_circ_native(hfss, "Setup1 : LastAdaptive", tag)
                res["outputs"].append(f"hfss_circ_{tag}.csv")
            except Exception as e:
                res["checks"]["circ_err"] = str(e)[:160]
                print(f"[!] {tag}: circular failed: {e}")

        pinfo = export_profile(hfss, None, tag, {"tag": tag, "k0a": v["k0a"], "mode": mode})
        res["outputs"].append(f"hfss_profile_{tag}.txt")
        res["peak_memory_MB"] = pinfo.get("peak_memory_MB")
        res["tetrahedra"] = pinfo.get("tetrahedra")
        hfss.save_project()
    finally:
        try:
            hfss.release_desktop(close_projects=False, close_desktop=True)
        except Exception:
            pass


# --------------------------------------------------- phase B: near-field map (big domain)
def phase_near(v, mode, tag, non_graphical, res):
    global AMP
    AMP = ANT_CAL if mode == "antenna" else 1.0
    hfss = _new_hfss(f"EX7_{tag}_near", non_graphical)
    try:
        hfss.modeler.model_units = "mm"
        build_body(hfss, v)
        _excite(hfss, v, mode, "x")
        _open_region(hfss, PAD_LAM)                             # big: reach +/-(a+2*lambda0)
        _setup(hfss, NEAR_PASSES)
        hfss.save_project()
        print(f"[i] {tag}: NEAR solve (big domain, passes={NEAR_PASSES}) ...")
        hfss.analyze_setup("Setup1", cores=CORES)
        n = export_near(hfss, "Setup1 : LastAdaptive", tag, v)
        res["outputs"].append(f"hfss_near_{tag}.csv ({n} pts)")
    finally:
        try:
            hfss.release_desktop(close_projects=False, close_desktop=True)
        except Exception:
            pass


# --------------------------------------------------------------------------- run one
def run_one(key, mode, non_graphical=True, do_circ=True, do_maps=True):
    v = next(x for x in CFG["variants"] if x["key"] == key)
    os.makedirs(OUT, exist_ok=True)
    tag = f"{key}_{mode}"
    res = {"tag": tag, "key": key, "mode": mode, "status": "?", "wall_s": None,
           "outputs": [], "checks": {}}
    print(f"\n[i] === {tag} : {v['title']} ({mode}) ===")
    t0 = time.perf_counter()
    try:
        phase_ff(v, mode, tag, non_graphical, do_circ, do_maps, res)
        res["status"] = "ok"
    except Exception as e:
        res["status"] = f"ff_failed: {type(e).__name__}: {str(e)[:180]}"
        print(f"[FAIL] {tag} FF: {res['status']}")
    if do_maps:
        try:
            phase_near(v, mode, tag, non_graphical, res)
        except Exception as e:
            res["checks"]["near_err"] = f"{type(e).__name__}: {str(e)[:180]}"
            print(f"[!] {tag}: near phase failed: {e}")
    res["wall_s"] = round(time.perf_counter() - t0, 1)
    return res


# --------------------------------------------------------------------------- driver
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Example 7 absolute fields (HFSS FEM)")
    ap.add_argument("selector", nargs="?", help="<key> | <body> | <key> <mode>")
    ap.add_argument("mode", nargs="?", choices=["diffraction", "antenna"])
    ap.add_argument("--graphical", action="store_true")
    ap.add_argument("--no-circ", action="store_true")
    ap.add_argument("--no-maps", action="store_true")
    ap.add_argument("--near-only", action="store_true",
                    help="re-run ONLY the near-field phase (re-solve near + Mag_E + vector Ex,Ey,Ez)")
    ap.add_argument("--modes", default="diffraction,antenna")
    a = ap.parse_args()

    keys = [v["key"] for v in CFG["variants"]]
    modes = [m for m in a.modes.split(",") if m]
    if a.selector in keys:
        sel = [a.selector]
        if a.mode:
            modes = [a.mode]
    elif a.selector in ("pec", "radome", "lens"):
        sel = [k for k in keys if k.startswith(a.selector + "_")]
    else:
        sel = keys

    results = []
    for k in sel:
        for m in modes:
            if a.near_only:                       # re-solve near only -> Mag_E + vector Ex,Ey,Ez
                v = next(x for x in CFG["variants"] if x["key"] == k)
                tag = f"{k}_{m}"
                res = {"tag": tag, "key": k, "mode": m, "outputs": []}
                print(f"\n[i] === NEAR-ONLY {tag} ===")
                try:
                    phase_near(v, m, tag, not a.graphical, res)
                    res["status"] = "ok"
                except Exception as e:
                    res["status"] = f"near_failed: {type(e).__name__}: {str(e)[:180]}"
                    print(f"[FAIL] {tag} near: {res['status']}")
                results.append(res)
                continue
            results.append(run_one(k, m, non_graphical=not a.graphical,
                                   do_circ=not a.no_circ, do_maps=not a.no_maps))
            with open(os.path.join(OUT, "ex7_run_summary.json"), "w", encoding="utf-8") as f:
                json.dump(results, f, indent=2, ensure_ascii=False)

    print("\n================ EX7 FEM SUMMARY ================")
    for r in results:
        print(f"{r['tag']:<28} {r.get('status'):<10} wall={r.get('wall_s')}s mem={r.get('peak_memory_MB')}MB "
              f"tets={r.get('tetrahedra')} out={len(r.get('outputs', []))}")
    print(f"\nOutputs in {os.path.relpath(OUT, HERE)}/  (summary: ex7_run_summary.json)")
