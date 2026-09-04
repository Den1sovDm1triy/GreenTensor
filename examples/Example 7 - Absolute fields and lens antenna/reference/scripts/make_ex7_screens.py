# SPDX-License-Identifier: MIT
"""Example 7 screenshots (graphical AEDT, native image exports). Model, cutaways, near-field
|E| heat maps, far-field 3D patterns, plane-wave setup, and SBR rays. Output: ex7_screens/.

Modes (curated, run what you need):
  geo   <key>          geometry isometric + cutaway (multilayer bodies)
  field <key>          diffraction: coarse solve -> near-field Mag_E map (x0z) + 3D RCS pattern
  ant   <key>          antenna: dipole + coarse solve -> 3D radiation pattern + near map
  sbr   <key>          SBR+: Visual Ray Trace (rays) + db(RCSTheta) cut
Uses native exporters (export_model_picture / export_field_jpg / export_report_to_jpg): the
AEDT window appears but images render to file (no foreground capture needed).
"""
import json
import os
import sys

from ansys.aedt.core import Hfss

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "ex7_screens")
os.makedirs(OUT, exist_ok=True)
CFG = json.load(open(os.path.join(HERE, "ex7_variants.json"), encoding="utf-8"))
AEDT_VERSION = os.environ.get("AEDT_VERSION", "2024.2")
CORES = int(os.environ.get("HFSS_CORES", "8"))
F_GHZ = CFG["frequency_Hz"] / 1e9
LAM0 = CFG["lambda0_mm"]
PASSES = int(os.environ.get("EX7_SHOT_PASSES", "4"))

COLORS = [[70, 130, 225], [60, 200, 120], [240, 180, 45], [210, 90, 60]]
PEC_COLOR = [160, 160, 170]
DIEL_COLOR = [90, 200, 210]


def _mat(hfss, mat):
    if mat["kind"] == "pec":
        return "pec"
    name = f"diel_e{mat['eps_re']:g}_t{mat['tand_e']:g}"
    if not hfss.materials.exists_material(name):
        m = hfss.materials.add_material(name)
        m.permittivity = mat["eps_re"]
        m.dielectric_loss_tangent = mat["tand_e"]
        m.update()
    return name


def build(hfss, v, colored=True):
    names, ci = [], 0
    for lyr in v["layers_hfss"]:
        mat = lyr["material"]
        if mat.get("kind") == "dielectric" and float(mat["eps_re"]) == 1.0 and float(mat["tand_e"]) == 0.0:
            continue
        mn = _mat(hfss, mat)
        obj = hfss.modeler.create_sphere([0, 0, 0], lyr["r_outer_mm"], name=f"L{lyr['index']}", material=mn)
        if lyr["r_inner_mm"] and lyr["r_inner_mm"] > 0:
            cut = hfss.modeler.create_sphere([0, 0, 0], lyr["r_inner_mm"], name=f"c{lyr['index']}")
            hfss.modeler.subtract(obj, cut, keep_originals=False)
            obj.material_name = mn
        if colored:
            obj.transparency = 0.0
            obj.color = PEC_COLOR if mat["kind"] == "pec" else (COLORS[ci % len(COLORS)] if mat["kind"] != "pec" else DIEL_COLOR)
            if mat["kind"] != "pec":
                ci += 1
        names.append(obj.name)
    return names


def _get(key):
    return next(x for x in CFG["variants"] if x["key"] == key)


def _new(design, sol="Modal"):
    return Hfss(design=design, solution_type=sol, version=AEDT_VERSION,
                non_graphical=False, new_desktop=True, close_on_exit=True)


def mode_geo(key):
    v = _get(key)
    hfss = _new(f"geo_{key}")
    try:
        hfss.modeler.model_units = "mm"
        names = build(hfss, v)
        hfss.modeler.fit_all()
        hfss.post.export_model_picture(os.path.join(OUT, f"{key}_model_iso.png"),
                                       orientation="isometric", show_axis=True, show_grid=False,
                                       show_ruler=False, width=1600, height=1000)
        if len(names) > 1:   # cutaway for multilayer (radome): remove y<0 half
            R = v["a_mm"] * 3 + 60
            box = hfss.modeler.create_box([-R, -R, -R], [2 * R, R, 2 * R], name="cutbox")
            hfss.modeler.subtract(names, box, keep_originals=False)
            hfss.modeler.fit_all()
            hfss.post.export_model_picture(os.path.join(OUT, f"{key}_cutaway.png"),
                                           orientation="isometric", show_axis=True, show_grid=False,
                                           show_ruler=False, width=1600, height=1000)
        print(f"[ok] geo {key}")
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


def _solve_diff(hfss, v):
    hfss.plane_wave(vector_format="Cartesian", origin=[0, 0, 0], polarization=[1, 0, 0],
                    propagation_vector=[0, 0, -1])
    pad_m = 1.25 * LAM0 / 1000.0
    f_op = 299792458.0 / (4.0 * pad_m)
    hfss.create_open_region(frequency=f"{f_op/1e9:.4f}GHz")
    hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=360, theta_step=2,
                                phi_start=0, phi_stop=90, phi_step=90, name="ff")
    s = hfss.create_setup("Setup1")
    s.props["Frequency"] = f"{F_GHZ}GHz"
    s.props["MaximumPasses"] = PASSES
    s.props["MaxDeltaS"] = 0.02
    s.update()
    hfss.analyze_setup("Setup1", cores=CORES)


def mode_field(key):
    v = _get(key)
    hfss = _new(f"fld_{key}")
    try:
        hfss.modeler.model_units = "mm"
        build(hfss, v, colored=False)
        _solve_diff(hfss, v)
        sol = "Setup1 : LastAdaptive"
        plot = hfss.post.create_fieldplot_cutplane(assignment="Global:XZ", quantity="Mag_E",
                                                   setup=sol, intrinsics={"Freq": f"{F_GHZ}GHz", "Phase": "0deg"},
                                                   plot_name=f"E_{key}")
        hfss.post.export_field_jpg(os.path.join(OUT, f"{key}_nearE.png"),
                                   getattr(plot, "name", f"E_{key}"), getattr(plot, "plot_folder", "") or "",
                                   orientation="front", width=1600, height=1000, display_wireframe=True,
                                   show_axis=True, show_grid=False, show_ruler=False)
        try:
            hfss.post.create_report(expressions="db(RCSTotal)", setup_sweep_name=sol,
                                    primary_sweep_variable="Theta", context="ff",
                                    report_category="Far Fields", plot_type="3D Polar Plot", plot_name="rcs3d")
            hfss.post.export_report_to_jpg(os.path.join(OUT, f"{key}_rcs3d.jpg"), "rcs3d", width=1400, height=1000)
        except Exception as e:
            print("   [rcs3d err]", e)
        print(f"[ok] field {key}")
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


def mode_ant(key):
    v = _get(key)
    l, w, g = LAM0 / 50.0, LAM0 / 400.0, LAM0 / 200.0
    a = v["a_mm"]
    hfss = _new(f"ant_{key}")
    try:
        hfss.modeler.model_units = "mm"
        build(hfss, v, colored=True)
        sh = hfss.modeler.create_rectangle("XZ", [-w / 2, 0.0, a + g], [l, w], name="dip")
        hfss.assign_current_source_to_sheet(sh.name, start_direction=[[0, 0, a + g], [0, 0, a + g + l]], name="I_z")
        pad_m = 1.25 * LAM0 / 1000.0
        hfss.create_open_region(frequency=f"{299792458.0/(4*pad_m)/1e9:.4f}GHz")
        hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=180, theta_step=2,
                                    phi_start=0, phi_stop=360, phi_step=10, name="ff")
        s = hfss.create_setup("Setup1")
        s.props["Frequency"] = f"{F_GHZ}GHz"
        s.props["MaximumPasses"] = PASSES
        s.update()
        hfss.modeler.fit_all()
        hfss.post.export_model_picture(os.path.join(OUT, f"{key}_antenna_model.png"),
                                       orientation="isometric", show_axis=True, show_grid=False,
                                       show_ruler=False, width=1600, height=1000)
        hfss.analyze_setup("Setup1", cores=CORES)
        sol = "Setup1 : LastAdaptive"
        try:
            hfss.post.create_report(expressions="db(DirTotal)", setup_sweep_name=sol,
                                    primary_sweep_variable="Theta", context="ff",
                                    report_category="Far Fields", plot_type="3D Polar Plot", plot_name="dir3d")
            hfss.post.export_report_to_jpg(os.path.join(OUT, f"{key}_antenna_pattern3d.jpg"), "dir3d",
                                           width=1400, height=1000)
        except Exception as e:
            print("   [pattern err]", e)
        print(f"[ok] ant {key}")
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


def mode_antcut(key):
    """Cross-section of the lens model WITH the inward Huygens source on the pole: full isometric
    then a half-cut (y<0 removed) revealing the ε=4 interior and the source sitting on the surface."""
    import run_ex7_antenna as ra
    v = _get(key)
    a = v["a_mm"]
    hfss = _new(f"antcut_{key}")
    try:
        hfss.modeler.model_units = "mm"
        names = build(hfss, v, colored=True)                 # lens sphere (coloured)
        z_surf = a
        if v["k0a"] >= 15:
            S = 1.6 * ra.DIP_L
            z_surf = ra.cut_facet(hfss, a, 1.5 * S, S / 2.0)
        ra.build_huygens(hfss, z_surf, ra.g_body(a), "H")    # inward Huygens element
        for nm, col in [("H_edip", [230, 40, 40]), ("H_mout", [245, 180, 40]), ("H_mbr", [245, 180, 40])]:
            try:
                o = hfss.modeler[nm]
                o.color = col
                o.transparency = 0.0
            except Exception:
                pass
        hfss.modeler.fit_all()
        hfss.post.export_model_picture(os.path.join(OUT, f"{key}_lens_src_iso.png"),
                                       orientation="isometric", show_axis=True, show_grid=False,
                                       show_ruler=False, width=1600, height=1000)
        # half-cut: remove y<0 of the LENS only (keep the source whole & visible on the cut face)
        R = a * 3 + 60
        box = hfss.modeler.create_box([-R, -R, -R], [2 * R, R, 2 * R], name="cutbox")
        hfss.modeler.subtract(names, box, keep_originals=False)
        hfss.modeler.fit_all()
        for orient in ("isometric", "front"):
            hfss.post.export_model_picture(os.path.join(OUT, f"{key}_lens_cutaway_{orient}.png"),
                                           orientation=orient, show_axis=True, show_grid=False,
                                           show_ruler=False, width=1600, height=1000)
        print(f"[ok] antcut {key}")
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


def mode_sbr(key):
    v = _get(key)
    hfss = _new(f"sbr_{key}", sol="SBR+")
    try:
        hfss.modeler.model_units = "mm"
        build(hfss, v, colored=True)
        hfss.plane_wave(vector_format="Cartesian", origin=[0, 0, 0], polarization=[1, 0, 0],
                        propagation_vector=[0, 0, -1])
        hfss.insert_infinite_sphere(definition="Theta-Phi", theta_start=0, theta_stop=360, theta_step=1,
                                    phi_start=0, phi_stop=0, phi_step=1, name="ff")
        st = hfss.create_setup("SBR")
        st.props["IsMonostaticRCS"] = False
        st.props["RadiationSetup"] = "ff"
        st.props["RayDensityPerWavelength"] = 10
        st.props["MaxNumberOfBounces"] = 5
        st.props["Frequency"] = f"{F_GHZ}GHz"
        st.update()
        hfss.modeler.fit_all()
        for mf, ith, tag in [("3GHz", 30, "obliq"), ("2GHz", 0, "head")]:
            try:
                vrt = hfss.post.create_sbr_plane_visual_ray_tracing(max_frequency=mf, ray_density=2,
                                                                    number_of_bounces=4, incident_theta=ith,
                                                                    incident_phi=0)
                hfss.modeler.fit_all()
                for orient in ["isometric", "right"]:
                    hfss.post.export_model_picture(os.path.join(OUT, f"{key}_rays_{tag}_{orient}.png"),
                                                   field_selections="all", orientation=orient,
                                                   show_axis=True, show_grid=False, show_ruler=False,
                                                   width=1600, height=1000)
                vrt.delete()
            except Exception as e:
                print(f"   [rays {tag} err]", e)
        try:
            hfss.analyze_setup("SBR", cores=CORES)
            sol = list(dict.fromkeys(hfss.existing_analysis_sweeps))[0]
            hfss.post.create_report(expressions="db(RCSTheta)", setup_sweep_name=sol,
                                    primary_sweep_variable="Theta", context="ff",
                                    variations={"Phi": "0deg"}, report_category="Far Fields",
                                    plot_type="Rectangular Plot", plot_name="rcscut")
            hfss.post.export_report_to_jpg(os.path.join(OUT, f"{key}_rcscut.jpg"), "rcscut", width=1400, height=800)
        except Exception as e:
            print("   [rcs err]", e)
        print(f"[ok] sbr {key}")
    finally:
        hfss.release_desktop(close_projects=True, close_desktop=True)


if __name__ == "__main__":
    mode, key = sys.argv[1], sys.argv[2]
    {"geo": mode_geo, "field": mode_field, "ant": mode_ant, "sbr": mode_sbr,
     "antcut": mode_antcut}[mode](key)
    print("Done:", mode, key)
