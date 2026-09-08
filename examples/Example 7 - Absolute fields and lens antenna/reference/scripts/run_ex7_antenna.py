# SPDX-License-Identifier: MIT
"""Example 7 ANTENNA mode, source = HUYGENS element (see ../MODELS.md).

Source at the +z pole = crossed TANGENTIAL electric dipole (a_theta = x, current sheet)
+ magnetic dipole (a_phi = y) realised as a small electric-current LOOP in the xz-plane
(PEC square washer with a gap fed by a current source). Own pattern = cardioid (1+cosθ)/2,
null at θ=180. Loop amplitude/phase are calibrated in FREE SPACE against that analytic
cardioid (NOT against GreenTensor), then FIXED for every body/size.

Diffraction mode is verified and untouched — this script only redoes the antenna branch.

Steps:
  python run_ex7_antenna.py --calib                 free-space acceptance test (§3)
  python run_ex7_antenna.py <key>|<body>|(all)      body antenna runs with the fixed source
  python run_ex7_antenna.py --standoff <key>        g/lambda0 convergence study (§4)
Flags: --graphical --no-circ --no-maps

Target free-space peak |E|(10 m) = I*k0*l*Z0/(2π r) = 0.75346 V/m (l=λ0/50, I=1 A) — this
is 2x the old vibrator (3.18) target 0.37673.
"""
import argparse
import json
import math
import os

import numpy as np

import run_ex7 as ex7                       # reuse verified helpers/config (loads CFG, utf-8 stdout)
from ansys.aedt.core import Hfss

OUT = ex7.OUT
CFG = ex7.CFG
F_GHZ = ex7.F_GHZ
LAM0 = ex7.LAM0
DIP_L = ex7.DIP_L
RE_SCALE = ex7.RE_SCALE
OBS_R = ex7.OBS_R
CORES = ex7.CORES
AEDT_VERSION = ex7.AEDT_VERSION
MAXPASSES = ex7.MAXPASSES
MAXDELTAS = ex7.MAXDELTAS

Z0 = CFG["antenna"]["Z0_ohm"]
K0 = CFG["k0_rad_m"]
E0_TARGET = 1.0 * K0 * (DIP_L / 1000.0) * Z0 / (2 * math.pi * OBS_R)   # 0.75346 V/m @10 m
W_SCALE = float(os.environ.get("EX7_WASHER_SCALE", "1.0"))            # magnetic-loop size scale
CALIB_JSON = os.path.join(OUT, "hfss_src_calib.json" if W_SCALE == 1.0
                          else f"hfss_src_calib_ws{W_SCALE:g}.json")
G_WORK = float(os.environ.get("EX7_G_LAM", str(1.0 / 400))) * LAM0     # working standoff (mm)


def g_body(a):
    """Working standoff on a body: max(fine g, 0.01*a). Large spheres need a meshable gap
    between the finite PEC source loop and the PEC surface (0.075 mm on a 143 mm domain hangs
    the mesher); 0.01*a is still a negligible fraction of the radius. Fine-g study is separate."""
    return max(G_WORK, 0.01 * a)


# --------------------------------------------------------------------- Huygens source
def build_huygens(hfss, a, g, prefix="H"):
    """Tangential x electric dipole (a_theta) + y magnetic dipole (a_phi) = Huygens element.
    Magnetic dipole = PEC square washer loop in the xz-plane (moment ‖ y), fed across a gap by
    a current source. A source-region mesh seed sharpens the small features. Returns
    (e_source_name, [magnetic source name])."""
    l, w = DIP_L, DIP_L / 8.0
    # Magnetic loop scale: EX7_WASHER_SCALE>1 enlarges the PEC loop (thicker bars, bigger gap,
    # more standoff) so it meshes on electrically-large domains (k15) without sliver wedges
    # against the sphere. The electric dipole (l, w) stays fixed (sets the absolute level); the
    # loop is re-calibrated in free space to the SAME cardioid, so the radiated field is identical.
    ws = float(os.environ.get("EX7_WASHER_SCALE", "1.0"))
    S, t, gp = 1.6 * l * ws, (l / 6.0) * ws, (l / 4.0) * ws
    # Whole source (co-located e-dipole + washer, half-height S/2) must clear the sphere:
    # put its BOTTOM at a+g so nothing dips below z=a into the body (fixes small-sphere crash).
    z0 = a + S / 2 + g
    # electric dipole: current sheet, current along x
    e_sh = hfss.modeler.create_rectangle("XY", [-l / 2, -w / 2, z0], [l, w], name=f"{prefix}_edip")
    e_il = [[-l / 2, 0.0, z0], [l / 2, 0.0, z0]]
    hfss.assign_current_source_to_sheet(e_sh.name, start_direction=e_il, name=f"{prefix}_Ie")
    # magnetic dipole: PEC washer loop (xz-plane) with a gap fed by a current source
    zc = z0
    outer = hfss.modeler.create_rectangle("XZ", [-S / 2, 0.0, zc - S / 2], [S, S], name=f"{prefix}_mout")
    inner = hfss.modeler.create_rectangle("XZ", [-(S / 2 - t), 0.0, zc - (S / 2 - t)],
                                          [S - 2 * t, S - 2 * t], name=f"{prefix}_min")
    hfss.modeler.subtract(outer, inner, keep_originals=False)
    gap = hfss.modeler.create_rectangle("XZ", [-gp / 2, 0.0, zc - S / 2], [t, gp], name=f"{prefix}_mgap")
    hfss.modeler.subtract(outer, gap, keep_originals=False)
    hfss.assign_perfecte_to_sheets([outer.name])
    br = hfss.modeler.create_rectangle("XZ", [-gp / 2, 0.0, zc - S / 2], [t, gp], name=f"{prefix}_mbr")
    m_il = [[-gp / 2, 0.0, zc - S / 2 + t / 2], [gp / 2, 0.0, zc - S / 2 + t / 2]]
    hfss.assign_current_source_to_sheet(br.name, start_direction=m_il, name=f"{prefix}_Im")
    return f"{prefix}_Ie", [f"{prefix}_Im"]


def cut_facet(hfss, a, half, depth):
    """Flatten a small facet at the +z pole (remove the curved cap over a square of half-width
    `half`, down to z=a-depth). The washer then sits on a FLAT surface -> the thin source/body
    gap has no curvature -> no sliver wedges. Returns the facet z-level."""
    z_flat = a - depth
    box = hfss.modeler.create_box([-half, -half, z_flat], [2 * half, 2 * half, a - z_flat + a],
                                  name="facetcut")
    bodies = [o for o in hfss.modeler.object_names if o.startswith("L")]
    hfss.modeler.subtract(bodies, box, keep_originals=False)
    print(f"[i] facet cut: half={half:.3g} mm, depth={depth:.3g} mm, z_flat={z_flat:.4g} mm")
    return z_flat


def _set(hfss, e_name, m_names, ie, im):
    """edit_sources: complex weight ie (electric) and im (each of the 4 loop segments)."""
    def ca(z):
        z = complex(z)
        return f"{abs(z)}A", f"{math.degrees(math.atan2(z.imag, z.real))}deg"
    d = {e_name: ca(ie)}
    for mn in m_names:
        d[mn] = ca(im)
    hfss.edit_sources(d)


def _open_region(hfss, pad_lam):
    pad_m = (pad_lam + 0.25) * LAM0 / 1000.0
    f_op = 299792458.0 / (4.0 * pad_m)
    hfss.create_open_region(frequency=f"{f_op/1e9:.4f}GHz")


# --------------------------------------------------------------------- free-space calib
def calibrate(non_graphical=True):
    os.makedirs(OUT, exist_ok=True)
    proj = "EX7A_calib"
    ex7._clean(proj)
    hfss = Hfss(project=os.path.join(OUT, proj + ".aedt"), design=proj, solution_type="Modal",
                version=AEDT_VERSION, non_graphical=non_graphical, new_desktop=True)
    try:
        hfss.modeler.model_units = "mm"
        e_name, m_names = build_huygens(hfss, 0.0, G_WORK, "H")   # source in vacuum at origin+g
        _open_region(hfss, 1.0)
        hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=180,
                                    theta_step=1, phi_start=0, phi_stop=90, phi_step=90, name="ff")
        st = hfss.create_setup("Setup1")
        st.props["Frequency"] = f"{F_GHZ}GHz"
        st.props["MaximumPasses"] = MAXPASSES
        st.props["MaxDeltaS"] = MAXDELTAS
        st.props["MinimumConvergedPasses"] = 2
        try:
            st.props["BasisOrder"] = 1
        except Exception:
            pass
        st.update()
        hfss.save_project()
        print("[calib] solving free-space Huygens source ...")
        hfss.analyze_setup("Setup1", cores=CORES)
        sol = "Setup1 : LastAdaptive"

        # electric-only and magnetic-only far fields (reuse one solve, weight via edit_sources)
        _set(hfss, e_name, m_names, 1.0, 0.0)
        th, ee0, ep_e0 = ex7._far_complex(hfss, sol, 0)
        _, ee90, ep_e90 = ex7._far_complex(hfss, sol, 90)
        _set(hfss, e_name, m_names, 0.0, 1.0)
        _, em0, ep_m0 = ex7._far_complex(hfss, sol, 0)
        _, em90, ep_m90 = ex7._far_complex(hfss, sol, 90)

        def at(arr, deg):
            return arr[int(np.argmin(np.abs(th - deg)))]
        # The source sits at the +z pole, so in BODY coordinates theta=0 points OUTWARD (away
        # from the body) and theta=180 points INTO the body. The feed must ILLUMINATE the body
        # the cardioid must point INWARD: null at theta=0
        # (outward), max at theta=180 (into body), shape (1-cos)/2. Least-squares null of the
        # OUTWARD radiation (theta=0) in BOTH planes.
        a, A = at(ee0, 0.0), at(em0, 0.0)
        b, B = at(ep_e90, 0.0), at(ep_m90, 0.0)
        c = -(np.conj(A) * a + np.conj(B) * b) / (abs(A) ** 2 + abs(B) ** 2)
        # total (relative), then scale so the INWARD peak |E|/10m = E0_TARGET
        tot_e0 = ee0 + c * em0
        peak = float(np.max(np.abs(tot_e0))) / OBS_R
        scale = E0_TARGET / peak if peak > 0 else 1.0

        # calibrated total fields
        te0 = scale * (ee0 + c * em0)
        te90_theta = scale * (ee90 + c * em90)      # rETheta at phi=90 (cross for antenna)
        tp90 = scale * (ep_e90 + c * ep_m90)        # rEPhi at phi=90 (co in H-plane)
        # metrics — cardioid now points INWARD: null at theta=0, max at theta=180
        pk = float(np.max(np.abs(te0)))
        null_db = 20 * math.log10(abs(at(te0, 0.0)) / pk + 1e-30)
        card = (1 - np.cos(np.radians(th))) / 2.0
        m = card > 0.05
        shp = 20 * np.log10(np.abs(te0)[m] / pk + 1e-30) - 20 * np.log10(card[m] + 1e-30)
        shape_med = float(np.median(np.abs(shp - np.median(shp))))
        level_db = 20 * math.log10((pk / OBS_R) / E0_TARGET)
        # phi0 vs phi90 co-pol agreement (Etheta@0 vs Ephi@90). Report over card>0.05 AND over
        # the significant lobes only (level >1% of peak, per the guide) — the deep null tail is
        # physically irrelevant and inflates dB ratios of two vanishing numbers.
        pkp = float(np.max(np.abs(tp90)))
        pp = 20 * np.log10(np.abs(te0) / pk + 1e-30) - 20 * np.log10(np.abs(tp90) / (pkp + 1e-30) + 1e-30)
        plane_med = float(np.median(np.abs(pp[m] - np.median(pp[m]))))
        sig = (np.abs(te0) / pk > 0.01) & (np.abs(tp90) / pkp > 0.01)
        plane_med_sig = float(np.median(np.abs(pp[sig] - np.median(pp[sig]))))

        ex7._write(os.path.join(OUT, "hfss_src_freespace.csv"),
                   "Theta,rETheta_phi0_V,rEPhi_phi90_V,rETheta_phi90_cross_V",
                   [th, np.abs(te0), np.abs(tp90), np.abs(te90_theta)])
        calib = {"scale_A": scale, "c_re": c.real, "c_im": c.imag, "c_mag": abs(c),
                 "c_phase_deg": math.degrees(math.atan2(c.imag, c.real)),
                 "E0_target_Vm": E0_TARGET, "peak_Vm_10m": pk / OBS_R,
                 "null_depth_dB": null_db, "shape_median_dB": shape_med,
                 "abs_level_err_dB": level_db, "plane_match_median_dB": plane_med,
                 "plane_match_sig_lobes_dB": plane_med_sig,
                 "g_work_mm": G_WORK, "dipole_len_mm": DIP_L}
        json.dump(calib, open(CALIB_JSON, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
        print(f"[calib] INWARD cardioid: c={abs(c):.4g}∠{calib['c_phase_deg']:.1f}°, scale={scale:.4g} A")
        print(f"[calib] null@0(outward)={null_db:.1f} dB (want ≤-30) | shape vs (1-cos)/2 "
              f"={shape_med:.2f} dB (≤0.3) | level err={level_db:+.2f} dB (≤0.5) | plane match="
              f"{plane_med:.2f} dB / sig-lobes {plane_med_sig:.2f} dB (≤0.3)")
        print(f"[calib] inward peak |E|(θ=180,10m)={pk/OBS_R:.4f} V/m (target {E0_TARGET:.4f})")
        ok = (null_db <= -30) and (shape_med <= 0.3) and (abs(level_db) <= 0.5) and (plane_med_sig <= 0.3)
        print(f"[calib] ACCEPTANCE (§3): {'PASS' if ok else 'FAIL — inspect before body runs'}")
        hfss.save_project()
        return calib
    finally:
        try:
            hfss.release_desktop(close_projects=False, close_desktop=True)
        except Exception:
            pass


def load_calib():
    if not os.path.isfile(CALIB_JSON):
        raise SystemExit("run --calib first (no hfss_src_calib.json)")
    d = json.load(open(CALIB_JSON, encoding="utf-8"))
    return d["scale_A"], complex(d["c_re"], d["c_im"])


# --------------------------------------------------------------------- body antenna run
def _setup_body(hfss, v):
    _open_region(hfss, ex7.PAD_LAM)
    hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=360, theta_step=1,
                                phi_start=0, phi_stop=270, phi_step=90, name="ff")
    hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=180, theta_step=3,
                                phi_start=0, phi_stop=360, phi_step=5, name="ffmap")
    st = hfss.create_setup("Setup1")
    st.props["Frequency"] = f"{F_GHZ}GHz"
    st.props["MaximumPasses"] = MAXPASSES
    st.props["MaxDeltaS"] = MAXDELTAS
    st.props["MinimumConvergedPasses"] = 2
    try:
        st.props["BasisOrder"] = 1
    except Exception:
        pass
    st.update()
    return st


def _circ_from_element1(hfss, sol, tag):
    """Circular = element1 + j*element2, element2 = element1 rotated 90° about z (sphere is
    z-symmetric): E_circ(θ,φ=0) = E1(θ,0) + j*E1(θ,270). rE_R/rE_L are circular projections."""
    th, et0, ep0 = ex7._far_complex(hfss, sol, 0)
    _, et270, ep270 = ex7._far_complex(hfss, sol, 270)
    ec_t = et0 + 1j * et270
    ec_p = ep0 + 1j * ep270
    s2 = math.sqrt(2.0)
    reR = np.abs((ec_t - 1j * ec_p) / s2)
    reL = np.abs((ec_t + 1j * ec_p) / s2)
    ex7._write(os.path.join(OUT, f"hfss_circ_{tag}.csv"), "Theta,rELHCP_V,rERHCP_V", [th, reL, reR])


def run_antenna(key, non_graphical=True, do_circ=True, do_maps=True):
    v = next(x for x in CFG["variants"] if x["key"] == key)
    scale, c = load_calib()
    tag = f"{key}_antenna"
    proj = f"EX7A_{key}"
    ex7._clean(proj)
    res = {"tag": tag, "key": key, "mode": "antenna", "source": "huygens", "status": "?", "outputs": []}
    print(f"\n[i] === {tag} : {v['title']} (Huygens) ===")
    hfss = Hfss(project=os.path.join(OUT, proj + ".aedt"), design=proj, solution_type="Modal",
                version=AEDT_VERSION, non_graphical=non_graphical, new_desktop=True)
    try:
        hfss.modeler.model_units = "mm"
        ex7.build_body(hfss, v)
        # Electrically-large sphere: flatten a facet at the pole so the flat PEC loop sits on a
        # flat surface (avoids the sliver mesh hang). Facet is tiny vs a -> optically negligible.
        z_surf = v["a_mm"]
        if v["k0a"] >= 15:
            S = 1.6 * DIP_L * W_SCALE
            z_surf = cut_facet(hfss, v["a_mm"], 1.5 * S, S / 2.0)
        e_name, m_names = build_huygens(hfss, z_surf, g_body(v["a_mm"]), "H")
        _setup_body(hfss, v)
        hfss.save_project()
        print(f"[i] {tag}: solving (cores={CORES}) ...")
        hfss.analyze_setup("Setup1", cores=CORES)
        _set(hfss, e_name, m_names, scale, scale * c)          # fixed calibrated Huygens weights
        sol = "Setup1 : LastAdaptive"

        th0, et0, _ = ex7._far_complex(hfss, sol, 0)
        _, et90, ep90 = ex7._far_complex(hfss, sol, 90)
        ex7._write(os.path.join(OUT, f"hfss_far_{tag}.csv"),
                   "Theta,rETheta_phi0_V,rEPhi_phi90_V,rETheta_phi90_cross_V",
                   [th0, np.abs(et0), np.abs(ep90), np.abs(et90)])
        res["outputs"].append(f"hfss_far_{tag}.csv")
        ref = v["reference_Vm"]["antenna"]
        pkE = float(np.max(np.abs(et0))) / OBS_R
        i90 = int(np.argmin(np.abs(th0 - 90.0)))
        e90 = float(np.abs(et0[i90])) / OBS_R
        ipk = int(np.argmax(np.abs(et0)))
        res["peak_Eplane_Vm"] = pkE
        res["peak_theta_deg"] = float(th0[ipk])
        res["E_theta90_Vm"] = e90
        _rf = ref['Eplane'] if ref['Eplane'] is not None else float('nan')
        print(f"[check] {tag}: peak |E|={pkE:.4g} V/m @θ={th0[ipk]:.0f}, broadside |E|(θ=90)={e90:.4g} V/m "
              f"(guide ref {_rf:.4g}); shape overlay is customer-side (§8)")

        if do_circ:
            try:
                _circ_from_element1(hfss, sol, tag)
                res["outputs"].append(f"hfss_circ_{tag}.csv")
            except Exception as e:
                res["circ_err"] = str(e)[:150]
                print("[!] circ:", e)
        if do_maps:
            for fn, lbl in ((lambda: ex7.export_near(hfss, sol, tag, v), "near"),
                            (lambda: ex7.export_farmap(hfss, sol, tag), "farmap")):
                try:
                    fn()
                    res["outputs"].append(f"hfss_{lbl}_{tag}.csv")
                except Exception as e:
                    res[f"{lbl}_err"] = str(e)[:150]
                    print(f"[!] {lbl}:", e)
        try:
            pinfo = ex7.export_profile(hfss, None, tag, {"tag": tag, "k0a": v["k0a"], "mode": "antenna"})
            res["outputs"].append(f"hfss_profile_{tag}.txt")
            res["peak_memory_MB"] = pinfo.get("peak_memory_MB")
            res["tetrahedra"] = pinfo.get("tetrahedra")
        except Exception as e:
            res["profile_err"] = str(e)[:150]
        res["status"] = "ok"
        hfss.save_project()
    except Exception as e:
        res["status"] = f"failed: {type(e).__name__}: {str(e)[:200]}"
        print(f"[FAIL] {tag}: {res['status']}")
    finally:
        try:
            hfss.release_desktop(close_projects=False, close_desktop=True)
        except Exception:
            pass
    return res


def run_standoff(key, non_graphical=True):
    """g/lambda0 convergence: peak level + shape median vs standoff (§4). One CSV per key."""
    v = next(x for x in CFG["variants"] if x["key"] == key)
    scale, c = load_calib()
    g_lams = [1 / 200.0, 1 / 400.0, 1 / 800.0]
    rows = []
    for gl in g_lams:
        g = gl * LAM0
        proj = f"EX7A_so_{key}"
        ex7._clean(proj)
        hfss = Hfss(project=os.path.join(OUT, proj + ".aedt"), design=proj, solution_type="Modal",
                    version=AEDT_VERSION, non_graphical=non_graphical, new_desktop=True)
        try:
            hfss.modeler.model_units = "mm"
            ex7.build_body(hfss, v)
            e_name, m_names = build_huygens(hfss, v["a_mm"], g, "H")
            _open_region(hfss, ex7.PAD_LAM)
            hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=360,
                                        theta_step=1, phi_start=0, phi_stop=90, phi_step=90, name="ff")
            st = hfss.create_setup("Setup1")
            st.props["Frequency"] = f"{F_GHZ}GHz"
            st.props["MaximumPasses"] = MAXPASSES
            st.props["MaxDeltaS"] = MAXDELTAS
            st.update()
            hfss.analyze_setup("Setup1", cores=CORES)
            _set(hfss, e_name, m_names, scale, scale * c)
            th, et0, _ = ex7._far_complex(hfss, "Setup1 : LastAdaptive", 0)
            pk = float(np.max(np.abs(et0))) / OBS_R
            rows.append((gl, pk))
            print(f"[standoff] {key} g/λ0=1/{1/gl:.0f}: peak |E|={pk:.4g} V/m")
        finally:
            try:
                hfss.release_desktop(close_projects=False, close_desktop=True)
            except Exception:
                pass
    ex7._write(os.path.join(OUT, f"hfss_standoff_{key}.csv"),
               "g_over_lambda0,peak_E_Vm", [[r[0] for r in rows], [r[1] for r in rows]])
    print(f"[OK] standoff -> hfss_standoff_{key}.csv")


def post_only(key, non_graphical=True):
    """Re-export circ/near/farmap/profile from an already-SOLVED antenna project (no re-solve).
    Recovers a run that crashed in post-processing after the (expensive) solve was saved."""
    v = next(x for x in CFG["variants"] if x["key"] == key)
    scale, c = load_calib()
    tag = f"{key}_antenna"
    proj = f"EX7A_{key}"
    hfss = Hfss(project=os.path.join(OUT, proj + ".aedt"), design=proj, solution_type="Modal",
                version=AEDT_VERSION, non_graphical=non_graphical, new_desktop=True)
    try:
        sol = "Setup1 : LastAdaptive"
        _set(hfss, "H_Ie", ["H_Im"], scale, scale * c)
        th0, et0, _ = ex7._far_complex(hfss, sol, 0)
        _, et90, ep90 = ex7._far_complex(hfss, sol, 90)
        ex7._write(os.path.join(OUT, f"hfss_far_{tag}.csv"),
                   "Theta,rETheta_phi0_V,rEPhi_phi90_V,rETheta_phi90_cross_V",
                   [th0, np.abs(et0), np.abs(ep90), np.abs(et90)])
        for fn, lbl in ((lambda: _circ_from_element1(hfss, sol, tag), "circ"),
                        (lambda: ex7.export_near(hfss, sol, tag, v), "near"),
                        (lambda: ex7.export_farmap(hfss, sol, tag), "farmap")):
            try:
                fn()
                print(f"[post] {tag}: {lbl} ok")
            except Exception as e:
                print(f"[post] {tag}: {lbl} FAILED: {e}")
        ex7.export_profile(hfss, None, tag, {"tag": tag, "k0a": v["k0a"], "mode": "antenna"})
        pkE = float(np.max(np.abs(et0))) / OBS_R
        print(f"[post] {tag}: peak |E|={pkE:.4g} V/m @θ={th0[int(np.argmax(np.abs(et0)))]:.0f}")
    finally:
        try:
            hfss.release_desktop(close_projects=False, close_desktop=True)
        except Exception:
            pass


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Example 7 antenna mode, Huygens element source")
    ap.add_argument("selector", nargs="?")
    ap.add_argument("--calib", action="store_true")
    ap.add_argument("--standoff", metavar="KEY")
    ap.add_argument("--post", metavar="KEY")
    ap.add_argument("--graphical", action="store_true")
    ap.add_argument("--no-circ", action="store_true")
    ap.add_argument("--no-maps", action="store_true")
    a = ap.parse_args()
    ng = not a.graphical

    if a.calib:
        calibrate(non_graphical=ng)
    elif a.post:
        post_only(a.post, non_graphical=ng)
    elif a.standoff:
        run_standoff(a.standoff, non_graphical=ng)
    else:
        keys = [v["key"] for v in CFG["variants"]]
        if a.selector in keys:
            sel = [a.selector]
        elif a.selector in ("pec", "radome", "lens"):
            sel = [k for k in keys if k.startswith(a.selector + "_")]
        elif a.selector == "small":     # all k0a=0.5,2 (k15 finite-source mesh hangs; documented)
            sel = [k for k in keys if not k.endswith("_k15")]
        elif a.selector == "large":
            sel = [k for k in keys if k.endswith("_k15")]
        else:
            sel = keys
        results = []
        for k in sel:
            results.append(run_antenna(k, non_graphical=ng, do_circ=not a.no_circ, do_maps=not a.no_maps))
            json.dump(results, open(os.path.join(OUT, "ex7_antenna_summary.json"), "w",
                                    encoding="utf-8"), indent=2, ensure_ascii=False)
        print("\n===== EX7 ANTENNA (Huygens) SUMMARY =====")
        for r in results:
            print(f"{r['tag']:<26} {r['status']:<10} peak={r.get('peak_Eplane_Vm')} "
                  f"mem={r.get('peak_memory_MB')}MB out={len(r['outputs'])}")
