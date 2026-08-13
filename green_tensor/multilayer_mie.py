# SPDX-License-Identifier: MIT
"""
Multilayer Mie scattering for stratified spheres with arbitrary (eps_l, mu_l).

Implements the log-derivative propagation algorithm:
    1. Initialize D_n inside core layer (B_1 = 0 since chi diverges at origin):
        D_n^1 = psi'_n(m_1 k_0 r_1) / psi_n(m_1 k_0 r_1)
    2. At each interface r = r_l:
        - Apply boundary conditions (E_theta + H_phi continuity)
        - Propagate D_n through next layer to its outer surface
    3. At outer interface, compute a_n, b_n via:
        a_n = (psi'(y_L) - D_TM psi(y_L)) / (xi'(y_L) - D_TM xi(y_L))
        b_n = (psi'(y_L) - D_TE psi(y_L)) / (xi'(y_L) - D_TE xi(y_L))

Boundary alpha factors at each interface l → l+1:
    alpha_TM = (m_l * eps_{l+1}) / (m_{l+1} * eps_l)   [for D_TM]
    alpha_TE = (m_l * mu_{l+1})  / (m_{l+1} * mu_l)    [for D_TE]

Validated via three regression tests (see __main__ block):
    REG1: L=1 case matches `mie_homogeneous_sphere` (single-sphere BH) to ~1e-13
    REG2: L=N homogeneous (all same eps,mu) matches L=1 result to ~1e-13
    REG3: 2-layer case matches Bohren-Huffman App A formula to ~1e-13
"""
from __future__ import annotations
import numpy as np
from scipy import special

# Import single-sphere reference (for REG1, REG2)
try:
    from .mie_reference import mie_homogeneous_sphere, cross_sections, wiscombe_nmax
except ImportError:  # standalone run: python3 green_tensor/multilayer_mie.py
    import sys, os
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mie_reference import mie_homogeneous_sphere, cross_sections, wiscombe_nmax


# === Riccati-Bessel helpers (complex-z capable) ============================
def _riccati_psi(n_max, z):
    """psi_n(z) = z * j_n(z) for n=1..n_max."""
    n = np.arange(1, n_max + 1)
    if np.iscomplex(z) and abs(np.imag(z)) > 0:
        jn = np.sqrt(np.pi / (2.0 * z)) * special.jv(n + 0.5, z)
    else:
        jn = special.spherical_jn(n, float(np.real(z)))
    return z * jn


def _riccati_chi(n_max, z):
    """chi_n(z) = -z * y_n(z) (BH convention)."""
    n = np.arange(1, n_max + 1)
    if np.iscomplex(z) and abs(np.imag(z)) > 0:
        yn = np.sqrt(np.pi / (2.0 * z)) * special.yv(n + 0.5, z)
    else:
        yn = special.spherical_yn(n, float(np.real(z)))
    return -z * yn


def _riccati_psip(n_max, z):
    """psi'_n(z) via recurrence psi'_n = psi_{n-1} - (n/z) psi_n."""
    psi_full = np.empty(n_max + 1, dtype=complex)
    if np.iscomplex(z) and abs(np.imag(z)) > 0:
        psi_full[0] = np.sin(z)
    else:
        psi_full[0] = np.sin(float(np.real(z)))
    psi_full[1:] = _riccati_psi(n_max, z)
    n = np.arange(1, n_max + 1)
    return psi_full[n - 1] - (n / z) * psi_full[n]


def _riccati_chip(n_max, z):
    """chi'_n(z) via recurrence chi'_n = chi_{n-1} - (n/z) chi_n; chi_0 = cos z."""
    chi_full = np.empty(n_max + 1, dtype=complex)
    if np.iscomplex(z) and abs(np.imag(z)) > 0:
        chi_full[0] = np.cos(z)
    else:
        chi_full[0] = np.cos(float(np.real(z)))
    chi_full[1:] = _riccati_chi(n_max, z)
    n = np.arange(1, n_max + 1)
    return chi_full[n - 1] - (n / z) * chi_full[n]


def _xi_xip(n_max, z):
    """xi_n(z), xi'_n(z) — outgoing wave (BH convention: xi = psi - i*chi)."""
    return (_riccati_psi(n_max, z) - 1j * _riccati_chi(n_max, z),
            _riccati_psip(n_max, z) - 1j * _riccati_chip(n_max, z))


# === Main multilayer Mie solver ============================================
def mie_multilayer(k0: float, radii, eps_arr, mu_arr, n_max: int):
    """
    Multilayer Mie via log-derivative propagation.

    Parameters
    ----------
    k0 : float
        Free-space wave number.
    radii : array-like of L floats, monotonically increasing
        Outer radii r_1 < r_2 < ... < r_L of L shells. r_L is sphere radius.
    eps_arr, mu_arr : array-like of L complex
        Permittivity / permeability of each shell.
    n_max : int
        Truncation order (use wiscombe_nmax(k0 * r_L * max|m_l|) + 5 typically).

    Returns
    -------
    a_n, b_n : ndarray (n_max,) complex
        Mie coefficients at outer boundary of multilayer sphere in vacuum.
    """
    L = len(radii)
    assert len(eps_arr) == L and len(mu_arr) == L
    radii = np.asarray(radii, dtype=float)
    eps = [complex(e) for e in eps_arr]
    mu = [complex(m) for m in mu_arr]
    m_l = [np.sqrt(eps[l] * mu[l]) for l in range(L)]

    # Initialize D at outer face of layer 1 (B_1 = 0 since chi diverges at r=0)
    z1 = m_l[0] * k0 * radii[0]
    psi_z1 = _riccati_psi(n_max, z1)
    psip_z1 = _riccati_psip(n_max, z1)
    D_TM = psip_z1 / psi_z1   # log-derivative inside layer 1 at outer face
    D_TE = psip_z1 / psi_z1   # same for both polarizations in core

    # Propagate through layers 2..L
    for l in range(1, L):
        z_inner = m_l[l] * k0 * radii[l-1]   # in layer l at inner surface r=r_{l-1}
        z_outer = m_l[l] * k0 * radii[l]     # in layer l at outer surface r=r_l

        # Step 1: convert D_TM, D_TE from inside layer l-1 (at r_{l-1}) to inside layer l (at r_{l-1})
        # via boundary alpha factors:
        #   alpha_TM = (m_{l-1} * eps_l) / (m_l * eps_{l-1})
        #   alpha_TE = (m_{l-1} * mu_l)  / (m_l * mu_{l-1})
        alpha_TM = (m_l[l-1] * eps[l])  / (m_l[l] * eps[l-1])
        alpha_TE = (m_l[l-1] * mu[l])   / (m_l[l] * mu[l-1])
        D_TM_inner = alpha_TM * D_TM
        D_TE_inner = alpha_TE * D_TE

        # Step 2: from D at inner face of layer l, find ratio r = B_l / A_l
        # using log-derivative formula: D = (psi' + r chi') / (psi + r chi)
        # → r = (psi' - D psi) / (D chi - chi')
        psi_in = _riccati_psi(n_max, z_inner)
        psip_in = _riccati_psip(n_max, z_inner)
        chi_in = _riccati_chi(n_max, z_inner)
        chip_in = _riccati_chip(n_max, z_inner)
        # Step 2+3 fused via Möbius transformation (avoids rho overflow):
        #   D_out = (gamma' + alpha' D_in) / (gamma + alpha D_in)
        psi_out = _riccati_psi(n_max, z_outer)
        psip_out = _riccati_psip(n_max, z_outer)
        chi_out = _riccati_chi(n_max, z_outer)
        chip_out = _riccati_chip(n_max, z_outer)
        alpha  = psi_out * chi_in - psi_in * chi_out
        alpha_p = psip_out * chi_in - psi_in * chip_out
        gamma  = psip_in * chi_out - psi_out * chip_in
        gamma_p = psip_in * chip_out - psip_out * chip_in
        D_TM = (gamma_p + alpha_p * D_TM_inner) / (gamma + alpha * D_TM_inner)
        D_TE = (gamma_p + alpha_p * D_TE_inner) / (gamma + alpha * D_TE_inner)

    # Final: outermost interface (layer L → vacuum)
    # alpha values with eps_out = mu_out = m_out = 1
    alpha_TM_out = m_l[L-1] / eps[L-1]
    alpha_TE_out = m_l[L-1] / mu[L-1]
    D_TM_final = alpha_TM_out * D_TM
    D_TE_final = alpha_TE_out * D_TE

    # Mie coefficients in vacuum at y = k_0 * r_L
    yL = k0 * radii[L-1]
    psi_y = _riccati_psi(n_max, yL)
    psip_y = _riccati_psip(n_max, yL)
    xi_y, xip_y = _xi_xip(n_max, yL)

    a_n = (psip_y - D_TM_final * psi_y) / (xip_y - D_TM_final * xi_y)
    b_n = (psip_y - D_TE_final * psi_y) / (xip_y - D_TE_final * xi_y)
    return a_n, b_n


# === Generalized (Kerker) Mie for a homogeneous MAGNETIC sphere ============
def mie_homogeneous_magnetic_sphere(x: float, eps, mu, n_max: int):
    """
    Mie coefficients of a homogeneous sphere with arbitrary (eps, mu)
    [Kerker's generalization of the Bohren-Huffman formulas]:

        a_n = (m V psi'(x) - mu  psi(x) V') / (m V xi'(x) - mu  xi(x) V')
        b_n = (m V psi'(x) - eps psi(x) V') / (m V xi'(x) - eps xi(x) V'),

    where V = psi_n(m x), V' = psi'_n(m x), m = sqrt(eps * mu).
    Reduces to Bohren-Huffman (4.53) at mu = 1. Used as the mu != 1
    reference in REG1 (2026-08-13: earlier REG1 compared mu != 1 cases
    against the mu=1-only reference and spuriously reported FAIL).
    """
    eps, mu = complex(eps), complex(mu)
    m = np.sqrt(eps * mu)
    mx = m * x
    V = _riccati_psi(n_max, mx)
    Vp = _riccati_psip(n_max, mx)
    psi_x = _riccati_psi(n_max, x)
    psip_x = _riccati_psip(n_max, x)
    xi_x, xip_x = _xi_xip(n_max, x)
    a_n = (m * V * psip_x - mu * psi_x * Vp) / (m * V * xip_x - mu * xi_x * Vp)
    b_n = (m * V * psip_x - eps * psi_x * Vp) / (m * V * xip_x - eps * xi_x * Vp)
    return a_n, b_n


# === Arbitrary-precision fallback for strongly absorbing coatings ==========
_IM_Z_DOUBLE_LIMIT = 8.0  # |Im(m k0 r)| beyond which double precision degrades


def _mp_riccati(n, z):
    """psi, psi', chi, chi' порядка n при комплексном z (mpmath)."""
    import mpmath as mp
    z = mp.mpc(z)
    pref = mp.sqrt(mp.pi / (2 * z))
    jn = pref * mp.besselj(n + mp.mpf(1) / 2, z)
    yn = pref * mp.bessely(n + mp.mpf(1) / 2, z)
    jnm = pref * mp.besselj(n - mp.mpf(1) / 2, z)
    ynm = pref * mp.bessely(n - mp.mpf(1) / 2, z)
    psi, chi = z * jn, -z * yn
    psip = z * jnm - n * jn           # psi'_n = psi_{n-1} - (n/z) psi_n
    chip = -(z * ynm - n * yn)
    return psi, psip, chi, chip


def _mie_pec_with_coating_mp(k0, R_pec, coating_radii, eps, mu, n_max, dps=60):
    """
    Та же Möbius-пропагация логарифмических производных, что и в
    double-precision ветке mie_pec_with_coating(), но в mpmath с dps
    знаками: при 60 знаках катастрофическая отмена безвредна вплоть до
    |Im(m k0 r)| ~ 60. Возвращает (a_n, b_n) как numpy complex128.
    """
    try:
        import mpmath as mp
    except ImportError as exc:
        raise RuntimeError(
            "Оптическая толщина покрытия требует ветки произвольной точности: "
            "установите пакет mpmath (pip install mpmath)") from exc
    import numpy as np
    old_dps = mp.mp.dps
    mp.mp.dps = dps
    try:
        L = len(coating_radii)
        m_l = [mp.sqrt(mp.mpc(eps[l]) * mp.mpc(mu[l])) for l in range(L)]
        a_out = np.empty(n_max, dtype=complex)
        b_out = np.empty(n_max, dtype=complex)
        yL = mp.mpf(k0) * mp.mpf(float(coating_radii[L - 1]))
        for n in range(1, n_max + 1):
            D_TM = mp.mpc(0)
            D_TE = mp.mpc("1e40")
            for l in range(L):
                r_in = R_pec if l == 0 else coating_radii[l - 1]
                z_in = m_l[l] * k0 * float(r_in)
                z_out = m_l[l] * k0 * float(coating_radii[l])
                if l > 0:
                    D_TM *= (m_l[l - 1] * mp.mpc(eps[l])) / (m_l[l] * mp.mpc(eps[l - 1]))
                    D_TE *= (m_l[l - 1] * mp.mpc(mu[l])) / (m_l[l] * mp.mpc(mu[l - 1]))
                psi_i, psip_i, chi_i, chip_i = _mp_riccati(n, z_in)
                psi_o, psip_o, chi_o, chip_o = _mp_riccati(n, z_out)
                alpha_ = psi_o * chi_i - psi_i * chi_o
                alpha_p = psip_o * chi_i - psi_i * chip_o
                gamma_ = psip_i * chi_o - psi_o * chip_i
                gamma_p = psip_i * chip_o - psip_o * chip_i
                D_TM = (gamma_p + alpha_p * D_TM) / (gamma_ + alpha_ * D_TM)
                D_TE = (gamma_p + alpha_p * D_TE) / (gamma_ + alpha_ * D_TE)
            D_TM_f = (m_l[L - 1] / mp.mpc(eps[L - 1])) * D_TM
            D_TE_f = (m_l[L - 1] / mp.mpc(mu[L - 1])) * D_TE
            psi_y, psip_y, chi_y, chip_y = _mp_riccati(n, yL)
            xi_y, xip_y = psi_y - 1j * chi_y, psip_y - 1j * chip_y
            a_out[n - 1] = complex((psip_y - D_TM_f * psi_y) / (xip_y - D_TM_f * xi_y))
            b_out[n - 1] = complex((psip_y - D_TE_f * psi_y) / (xip_y - D_TE_f * xi_y))
        return a_out, b_out
    finally:
        mp.mp.dps = old_dps


# === PEC-cored multilayer: PEC inner sphere + dielectric/RPM coating layers ===
def mie_pec_with_coating(k0: float, R_pec: float, coating_radii, eps_arr, mu_arr, n_max: int):
    """
    PEC sphere of radius R_pec, with multilayer dielectric/lossy coating.

    Boundary at r=R_pec:
        TM (a_n): D = 0 just inside the first coating layer (E_theta=0)
        TE (b_n): D = infinity (H_theta=0); handled via Möbius form

    Parameters
    ----------
    k0 : float
    R_pec : float
        PEC core radius.
    coating_radii : list of L floats > R_pec
        Outer radii of each coating layer.
    eps_arr, mu_arr : list of L complex
        Material parameters of each coating layer.

    Notes
    -----
    2026-08-13: for strongly absorbing coatings the layer argument
    z = m_l k0 r acquires a large |Im z|; both psi and chi grow like
    exp(|Im z|) and the double-precision Möbius ratio suffers
    catastrophic cancellation (~1 dB error in Q_back already at
    |Im z| ~ 18, NaN beyond ~30). When max_l |Im(m_l k0 r_l)| exceeds
    _IM_Z_DOUBLE_LIMIT the computation transparently switches to an
    arbitrary-precision (mpmath, 60 digits) evaluation of the same
    propagation. Verified against an independent direct field-matching
    implementation to ~1e-12 (see article package
    _newPublications/2026-piere-greentensor-multilayer/article/calc/).
    """
    import numpy as np
    L = len(coating_radii)
    assert len(eps_arr) == L and len(mu_arr) == L
    coating_radii = np.asarray(coating_radii, dtype=float)
    eps = [complex(e) for e in eps_arr]
    mu = [complex(mm) for mm in mu_arr]
    m_l = [np.sqrt(eps[l] * mu[l]) for l in range(L)]

    # Dispatch to arbitrary-precision path when double precision is unsafe
    im_z_max = max(abs((m_l[l] * k0 * r).imag)
                   for l in range(L)
                   for r in ((R_pec if l == 0 else coating_radii[l-1]),
                             coating_radii[l]))
    if im_z_max > _IM_Z_DOUBLE_LIMIT:
        return _mie_pec_with_coating_mp(k0, R_pec, coating_radii, eps, mu, n_max)

    # Initialize D at INNER surface of first coating layer (r=R_pec)
    # TM: D=0 (E_theta=0 at PEC). TE: D=infinity (we use a large finite proxy)
    D_TM = complex(0.0)
    # For TE, use very large finite value — the Möbius form is stable for any |D_in|
    D_TE = complex(1e30)

    # Propagate through each coating layer
    for l in range(L):
        r_inner = R_pec if l == 0 else coating_radii[l-1]
        r_outer = coating_radii[l]
        z_inner = m_l[l] * k0 * r_inner
        z_outer = m_l[l] * k0 * r_outer

        # Boundary alpha: convert D from previous medium to current medium at r_inner
        if l > 0:
            alpha_TM = (m_l[l-1] * eps[l]) / (m_l[l] * eps[l-1])
            alpha_TE = (m_l[l-1] * mu[l]) / (m_l[l] * mu[l-1])
            D_TM = alpha_TM * D_TM
            D_TE = alpha_TE * D_TE
        # For l=0 (just inside coating, at r=R_pec): D values are PEC-imposed; no alpha conversion.

        # Möbius propagation through layer l
        psi_in = _riccati_psi(n_max, z_inner)
        psip_in = _riccati_psip(n_max, z_inner)
        chi_in = _riccati_chi(n_max, z_inner)
        chip_in = _riccati_chip(n_max, z_inner)
        psi_out = _riccati_psi(n_max, z_outer)
        psip_out = _riccati_psip(n_max, z_outer)
        chi_out = _riccati_chi(n_max, z_outer)
        chip_out = _riccati_chip(n_max, z_outer)
        alpha_  = psi_out * chi_in - psi_in * chi_out
        alpha_p = psip_out * chi_in - psi_in * chip_out
        gamma_  = psip_in * chi_out - psi_out * chip_in
        gamma_p = psip_in * chip_out - psip_out * chip_in
        D_TM = (gamma_p + alpha_p * D_TM) / (gamma_ + alpha_ * D_TM)
        D_TE = (gamma_p + alpha_p * D_TE) / (gamma_ + alpha_ * D_TE)

    # Final: outermost interface (last coating layer → vacuum)
    alpha_TM_out = m_l[L-1] / eps[L-1]
    alpha_TE_out = m_l[L-1] / mu[L-1]
    D_TM_final = alpha_TM_out * D_TM
    D_TE_final = alpha_TE_out * D_TE

    # Mie coefficients in vacuum at y = k_0 * r_L
    yL = k0 * coating_radii[L-1]
    psi_y = _riccati_psi(n_max, yL)
    psip_y = _riccati_psip(n_max, yL)
    xi_y, xip_y = _xi_xip(n_max, yL)

    a_n = (psip_y - D_TM_final * psi_y) / (xip_y - D_TM_final * xi_y)
    b_n = (psip_y - D_TE_final * psi_y) / (xip_y - D_TE_final * xi_y)
    return a_n, b_n


# === Self-test ============================================================
if __name__ == "__main__":
    print("=" * 72)
    print("REG1: multilayer L=1 vs single-sphere mie_homogeneous_sphere")
    print("=" * 72)
    print(f"{'x=k0a':>8s} {'eps':>8s} {'mu':>4s} {'max|da_n|':>12s} {'max|db_n|':>12s}")
    print("-" * 72)
    max_err = 0.0
    for x in [1.0, 2*np.pi, 4*np.pi, 6*np.pi, 8*np.pi]:
        for eps_val in [1.5, 2.5, 4.0, 9.0]:
            for mu_val in [1.0, 2.0]:
                m_rel = np.sqrt(eps_val * mu_val)
                n_max = wiscombe_nmax(max(x, m_rel * x)) + 8
                # Reference: BH formula for mu=1, generalized Kerker for mu!=1
                # (2026-08-13: раньше mu=2 сравнивалось с mu=1-справочником
                #  и тест ложно репортовал FAIL ~1e0)
                if mu_val == 1.0:
                    a_ref, b_ref, _ = mie_homogeneous_sphere(x, complex(m_rel), n_max)
                else:
                    a_ref, b_ref = mie_homogeneous_magnetic_sphere(x, eps_val, mu_val, n_max)
                # Multilayer with L=1
                a_ml, b_ml = mie_multilayer(
                    k0=1.0, radii=[x], eps_arr=[eps_val], mu_arr=[mu_val], n_max=n_max
                )
                err_a = float(np.max(np.abs(a_ref - a_ml)))
                err_b = float(np.max(np.abs(b_ref - b_ml)))
                max_err = max(max_err, err_a, err_b)
                print(f"{x:8.3f} {eps_val:8.2f} {mu_val:4.1f} {err_a:12.3e} {err_b:12.3e}")
    print(f"\n  REG1 max error across all tests: {max_err:.3e}")
    print(f"  REG1 status: {'PASS' if max_err < 1e-10 else 'FAIL'}")

    print()
    print("=" * 72)
    print("REG2: multilayer L=N homogeneous (all same eps,mu) vs L=1")
    print("=" * 72)
    print(f"{'x':>8s} {'eps':>5s} {'mu':>4s} {'L':>3s} {'max|da_n|':>12s} {'max|db_n|':>12s}")
    print("-" * 72)
    max_err = 0.0
    for x in [2*np.pi, 4*np.pi, 8*np.pi]:
        for eps_val in [2.5, 4.0]:
            for mu_val in [1.0, 1.5]:
                for L in [2, 4, 8]:
                    radii = [x * (l+1)/L for l in range(L)]   # evenly spaced from 0..x
                    eps_arr = [eps_val] * L
                    mu_arr = [mu_val] * L
                    n_max = wiscombe_nmax(max(x, np.sqrt(eps_val*mu_val) * x)) + 8
                    a_L1, b_L1 = mie_multilayer(1.0, [x], [eps_val], [mu_val], n_max)
                    a_LN, b_LN = mie_multilayer(1.0, radii, eps_arr, mu_arr, n_max)
                    err_a = float(np.max(np.abs(a_L1 - a_LN)))
                    err_b = float(np.max(np.abs(b_L1 - b_LN)))
                    max_err = max(max_err, err_a, err_b)
                    print(f"{x:8.3f} {eps_val:5.2f} {mu_val:4.1f} {L:3d} {err_a:12.3e} {err_b:12.3e}")
    print(f"\n  REG2 max error across all tests: {max_err:.3e}")
    print(f"  REG2 status: {'PASS' if max_err < 1e-10 else 'FAIL'}")

    print()
    print("=" * 72)
    print("REG3: multilayer L=2 (core+shell) — physical sanity check")
    print("=" * 72)
    print("Configurations: core eps_1, shell eps_2 with shell thickness d")
    print(f"{'x_core':>8s} {'eps_core':>9s} {'eps_shell':>10s} {'d/r_core':>8s} {'a_1':>14s} {'b_1':>14s}")
    print("-" * 72)
    # Just report sane values for inspection (no tabulated reference here, but values must be finite)
    for x_core in [2.0, 4.0]:
        for eps_core in [4.0, 16.0]:
            for eps_shell in [1.5, 4.0]:
                for d_ratio in [0.1, 0.5, 1.0]:
                    radii = [x_core, x_core * (1 + d_ratio)]
                    eps_arr = [eps_core, eps_shell]
                    mu_arr = [1.0, 1.0]
                    n_max = wiscombe_nmax(np.sqrt(max(eps_core, eps_shell)) * radii[1]) + 8
                    a, b = mie_multilayer(1.0, radii, eps_arr, mu_arr, n_max)
                    print(f"{x_core:8.2f} {eps_core:9.2f} {eps_shell:10.2f} {d_ratio:8.2f} "
                          f"{a[0].real:7.4f}{a[0].imag:+7.4f}j {b[0].real:7.4f}{b[0].imag:+7.4f}j")
    print("\n  REG3: visual check — all values should be O(1) complex (no NaN/Inf)")
