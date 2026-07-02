# SPDX-License-Identifier: MIT

"""Единый раннер тестов GreenTensor без зависимости от pytest.

    python3 tests/run_all.py

Если установлен pytest, эквивалентно: pytest tests/.
"""
from __future__ import annotations

import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import analytic_mie  # noqa: E402
import test_sphere_mie  # noqa: E402
import test_sphere_multilayer as ml  # noqa: E402
import test_sphere_complex as cx  # noqa: E402
import test_sphere_magnetic as smg  # noqa: E402
import test_tmatrix as tmx  # noqa: E402
import test_vswf as vsw  # noqa: E402
import test_gmm as gmt  # noqa: E402
import test_cylinder as cyt  # noqa: E402
import test_solvers as slv  # noqa: E402
import test_cylinder_layered as cyl_l  # noqa: E402
import test_hfss_reference as hfss  # noqa: E402


def main() -> int:
    suite = [
        ("analytic_mie self-check (Rayleigh)", lambda: _run_module_main(analytic_mie)),
        ("sphere Q_sca vs Mie", test_sphere_mie.test_sphere_qsca_matches_analytic_mie),
        ("subdivision invariance", ml.test_subdivision_invariance),
        ("vacuum shell = core", ml.test_vacuum_shell_equals_core),
        ("backscatter vs Mie", ml.test_backscatter_matches_mie),
        ("absorption Q_sca/Q_abs vs Mie", cx.test_absorption_matches_mie),
        ("absorption sign (Im eps)", cx.test_absorption_sign),
        ("real metal vs Mie", cx.test_real_metal_matches_mie),
        ("PEC limit finite & correct", cx.test_pec_limit_is_finite_and_correct),
        ("metal shell shields -> PEC", cx.test_metal_shell_shields_to_pec),
        ("mag sphere formula reduces to nonmag", smg.test_eps_mu_formula_reduces_to_nonmagnetic_mie),
        ("mag sphere coeffs vs eps/mu Mie", smg.test_homogeneous_magnetodielectric_coefficients_match_closed_mie),
        ("mag sphere Qsca/Qext vs eps/mu Mie", smg.test_homogeneous_magnetodielectric_cross_sections_match_closed_mie),
        ("mag sphere distinct E/M response", smg.test_mu_only_sphere_has_distinct_electric_and_magnetic_response),
        ("mag sphere subdivision invariance", smg.test_magnetodielectric_subdivision_invariance),
        ("mag core metal shell -> PEC", smg.test_magnetic_core_metal_shell_shields_to_pec),
        ("T-matrix sphere vs analytic", tmx.test_sphere_tmatrix_cross_sections_vs_analytic),
        ("T-matrix optical theorem", tmx.test_from_ab_optical_theorem),
        ("T-matrix adapter consistency", tmx.test_adapter_matches_direct_construction),
        ("T-matrix param validation", tmx.test_tmatrix_validation),
        ("vswf Wigner-3j", vsw.test_wigner3j_known),
        ("vswf Gaunt vs quadrature", vsw.test_gaunt_vs_quadrature),
        ("vswf Y_nm orthonormal", vsw.test_ynm_orthonormal),
        ("vswf scalar addition theorem", vsw.test_scalar_addition_theorem),
        ("vswf zero-translation identity", vsw.test_zero_translation_identity),
        ("vswf vector M,N curl relations", vsw.test_vsw_curl_relations),
        ("vswf vector translation A,B", vsw.test_vector_translation),
        ("vswf closed-form translation (Cruzan)", vsw.test_translation_closed),
        ("gmm decoupling limit", gmt.test_decoupling_limit),
        ("gmm coupling is active", gmt.test_coupling_is_active),
        ("gmm solver residual", gmt.test_solver_residual),
        ("gmm single-sphere cross-sections", gmt.test_single_sphere_cross_sections),
        ("gmm cluster energy conservation", gmt.test_energy_conservation_lossless_cluster),
        ("gmm dense packing (closed Cruzan)", gmt.test_dense_packing_energy),
        ("gmm full-T == diagonal (regression)", gmt.test_fullT_matches_diagonal),
        ("gmm t_matrix path == t_vector path", gmt.test_t_matrix_path_matches_t_vector),
        ("gmm axial incidence finite (NaN fix)", gmt.test_axial_incidence_finite),
        ("gmm PEC boundary condition (phase)", gmt.test_pec_boundary_condition),
        ("gmm asymmetric cluster energy", gmt.test_asymmetric_cluster_energy),
        ("gmm overlap guard (ValueError)", gmt.test_overlap_guard),
        ("cylinder energy conservation", cyt.test_energy_conservation_lossless),
        ("cylinder absorption sign", cyt.test_absorption_sign),
        ("cylinder small-x behavior", cyt.test_small_x_positive_decreasing),
        ("cylinder finite not-impl", cyt.test_finite_not_implemented),
        ("solvers: SphereSolver vs canonical sphere", slv.test_sphere_solver_matches_core),
        ("solvers: sphere T-matrix consistency", slv.test_sphere_tmatrix_consistency),
        ("solvers: sphere eps/mu public API", slv.test_sphere_solver_magnetodielectric_matches_closed_mie),
        ("solvers: as_scatterer in Cluster", slv.test_sphere_as_scatterer_in_cluster),
        ("solvers: Cluster.solve shapes", slv.test_cluster_solve_shapes),
        ("solvers: CylinderSolver vs module", slv.test_cylinder_solver_matches_module),
        ("solvers: NotImplementedError honesty", slv.test_notimplemented_honesty),
        ("solvers: sphere antenna problem", slv.test_sphere_antenna_problem),
        ("sphere_mag: eps/mu formula -> nonmagnetic", smg.test_eps_mu_formula_reduces_to_nonmagnetic_mie),
        ("sphere_mag: coefficients vs closed Mie", smg.test_homogeneous_magnetodielectric_coefficients_match_closed_mie),
        ("sphere_mag: cross sections vs closed Mie", smg.test_homogeneous_magnetodielectric_cross_sections_match_closed_mie),
        ("sphere_mag: mu-only response", smg.test_mu_only_sphere_has_distinct_electric_and_magnetic_response),
        ("sphere_mag: dual eps=mu zero backscatter", smg.test_dual_sphere_zero_backscatter_and_equal_coefficients),
        ("sphere_mag: lossless eps,mu energy", smg.test_lossless_magnetodielectric_energy_conservation),
        ("sphere_mag: identical magnetic layers", smg.test_magnetodielectric_subdivision_invariance),
        ("sphere_mag: magnetic core + metal shell -> PEC", smg.test_magnetic_core_metal_shell_shields_to_pec),
        ("layered cyl: N=1 vs B&H arbiter", cyl_l.test_layered_n1_vs_arbiter),
        ("layered cyl: subdivision invariance", cyl_l.test_subdivision_invariance),
        ("layered cyl: energy conservation", cyl_l.test_energy_conservation_lossless),
        ("layered cyl: coated absorption", cyl_l.test_coated_absorption),
        ("layered cyl: solver wrappers match", cyl_l.test_solver_wrapper_matches),
        ("layered cyl: 1-layer vs homogeneous", cyl_l.test_solver_single_layer_matches_homogeneous),
        ("oblique cyl: N=1 vs Wait arbiter", cyl_l.test_oblique_n1_vs_arbiter),
        ("oblique cyl: reduces to normal", cyl_l.test_oblique_reduces_to_normal),
        ("oblique cyl: subdivision invariance", cyl_l.test_oblique_subdivision_invariance),
        ("oblique cyl: energy conservation", cyl_l.test_oblique_energy_conservation),
        ("oblique cyl: cross-sec reduce + wrappers", cyl_l.test_oblique_cross_sections_reduce_and_wrapper),
        ("sphere pattern vs Ansys HFSS (FEM)", hfss.test_pattern_vs_hfss_reference),
    ]
    failures = 0
    for name, fn in suite:
        try:
            fn()
            print(f"  ok   {name}")
        except Exception:  # noqa: BLE001
            failures += 1
            print(f"  FAIL {name}")
            traceback.print_exc()
    print(f"\n{'ВСЕ ТЕСТЫ ПРОЙДЕНЫ' if not failures else f'ПРОВАЛОВ: {failures}'}")
    return 1 if failures else 0


def _run_module_main(mod):
    # analytic_mie выполняет самопроверку в своём __main__; повторим ассерты тут
    m = 1.5
    for x in (0.01, 0.02, 0.05):
        exact = mod.q_sca(m, x, nmax=4)
        ray = mod.q_sca_rayleigh(m, x)
        assert abs(exact - ray) / ray < 0.03


if __name__ == "__main__":
    sys.exit(main())
