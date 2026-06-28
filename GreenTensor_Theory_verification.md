<!-- SPDX-License-Identifier: MIT -->

# Верификация GreenTensor_Theory.tex по публикациям

Каждый раздел независимо сверен с литературой (агент-скептик, web-поиск первоисточников).

> **Область охвата (v0.5.0).** В библиотеке и в `GreenTensor_Theory.tex` оставлены только
> точные аналитические ТФГ-разделы: слоистая сфера, бесконечный/слоистый/косой цилиндр и
> GMM-суперпозиция кластера сфер. Разделы «Конфокальный сфероид», «Триаксиальный эллипсоид»,
> «Конечный цилиндр», «Конечный конус», «EBCM» и «Разложение в сферы» **удалены** из теории
> в 0.5.0; соответствующие строки ниже — исторический аудит прежней редакции.

## Сводка

| Раздел | Вердикт | Подтв. | Испр. | Не пров. |
|---|---|---:|---:|---:|
| Постановка и архитектура | corrected | 7 | 4 | 0 |
| Слоистая сфера (ядро) | corrected | 13 | 2 | 0 |
| Конфокальный сфероид | corrected | 15 | 1 | 0 |
| Триаксиальный эллипсоид | corrected | 23 | 3 | 0 |
| Конечный цилиндр (полу-аналитика) | corrected | 13 | 1 | 0 |
| Конечный конус (полу-аналитика) | corrected | 11 | 2 | 1 |
| Суперпозиция T-матриц (GMM) | corrected | 14 | 3 | 1 |
| Разложение и верификация | corrected | 12 | 1 | 0 |

## Детали по разделам

### Постановка и архитектура
**Исправлено:**
- Sphere T-matrix diagonal T = diag(-a_n,-b_n), Mie coefficients; minus sign per S=I+2T convention  
  _источник:_ Waterman (1971), Bohren & Huffman (1983): sphere T diagonal = Mie coefficients; added explicit S=I+2T note to justify the minus-sign convention
- Lippmann-Schwinger volume integral: E = E_inc + int G0 . [k^2(r')-k_ext^2] E d^3r' with curl curl G0 - k_ext^2 G0 = I delta  
  _источник:_ Standard EM volume integral equation (Tsang/Kong, Mishchenko 2002). Rewrote source term as exact contrast [k^2-k_ext^2] valid for magnetodielectric media; original omega^2 mu_0 Delta-eps form retained as the non-magnetic reduction (it was correct only for mu=mu_0 despite prose claiming general magnetodielectric)
- Eisenhart/Stackel: scalar Helmholtz separates in exactly 11 orthogonal coordinate systems (confocal quadrics and degenerations); sphere/prolate/oblate spheroidal/confocal ellipsoidal among them  
  _источник:_ Eisenhart (1934); Weinacht-Eisenhart theorem. Refined wording to 'конфокальных квадрик и их вырождений'; removed slightly inaccurate 'R-разделяется' (Helmholtz with k!=0 is simple-separable, not R-separable) - R-separation applies to special cases only
- Spheroid confocal condition a_j^2 - b_j^2 = (d/2)^2 (constant semifocal distance per shell), aspect ratio varies shell-to-shell  
  _источник:_ Farafonov & Voshchinnikov (2012); Farafonov et al. (1996). Confirmed a^2-b^2=const=semifocal-distance^2. Clarified that d is full inter-focal distance and d/2 the semifocal, removing ambiguity
**Открытые оговорки:**
- The Lippmann-Schwinger term was rewritten to the exact contrast form [k^2(r')-k_ext^2]E, which is rigorous for magnetodielectric media; the original omega^2 mu_0 Delta-epsilon form is now presented as the explicit non-magnetic (mu=mu_0) reduction. Strictly, for a magnetic scatterer the full volume integral also contains a magnetic-contrast term involving curl operations on G0; this is omitted because the section uses the integral equation only as conceptual motivation and the analytic solvers do not solve it directly. If full magnetodielectric rigor is desired in this passage, add the magnetic-current/contrast term explicitly.
- eq:adenkerker uses scalar epsilon-only impedance jumps sqrt(eps_{h+1}/eps_h); this matches the non-magnetic core in mie_core.py. The code and prose elsewhere allow complex mu per layer (sigma = sqrt(eps*mu)). For a fully magnetic layered sphere the TM/TE jump factors and Bessel arguments carry mu; the displayed recursion is the non-magnetic specialization. Consider noting mu-dependence for full generality.
- The claim that 10 of the 11 separable systems are degenerations of confocal ellipsoidal (with paraboloidal as exception) is correct but not stated; current wording 'конфокальных квадрик и их вырождений' is sufficient and supported by Eisenhart (1934).
- Citation flammer1957, vanburen2020, dlmf are not used in this section (they belong to later spheroidal-detail sections); this is acceptable for a postановка/architecture section but flagged for the shared-bibliography editor to confirm they are cited elsewhere in the document.

### Слоистая сфера (ядро)
**Исправлено:**
- Spherical Bessel derivative recurrence psi_n'(x) = n/(2n+1) psi_{n-1} - (n+1)/(2n+1) psi_{n+1} + psi_n/x  
  _источник:_ DLMF 10.51.2: j_n'(x) = n/(2n+1) j_{n-1} - (n+1)/(2n+1) j_{n+1}, times x plus psi_n/x. Original LaTeX had a garbled psi_k^flat definition (x j_k/x . x = x j_k, nonsense); replaced with clean recurrence matching code _psi_chi.
- tau_n identity via Legendre P_n^0, P_n^2  
  _источник:_ scipy.special.lpmv docs confirm Condon-Shortley phase (-1)^m IS included. With CS phase, dP_n^1/dtheta = (1/2)[n(n+1)P_n^0 - P_n^2], so code's tay = (1/2)[P_n^2 - n(n+1)P_n^0] = -tau_n. Original LaTeX claimed exact equivalence to tau_n; corrected to state tay = -tau_n with explicit CS-phase note (sign cancels via modulus).
**Открытые оговорки:**
- Identification M_n,N_n <-> Mie a_n,b_n is only up to a phase/sign convention fixed by the code's xi_n scaling and the conjugation (np.conj on lines 168). The original claim 'M_n=-a_n, N_n=-b_n' was internally contradictory with Q_ext = (2/x^2) sum (2n+1) Re(M_n+N_n) > 0 (a global -1 would flip the sign of extinction). I softened the statement to 'equal up to phase convention, chosen so Q_ext>0, verified against analytic Mie' rather than asserting an exact sign. If a precise sign relation is desired, it should be derived from the code's actual conjugation+scaling, not stated as -a_n.
- Notational mismatch (not a physics error): the document defines chi_n = -x y_n (Bohren-Huffman), but the code stores the singular solution as +x y_n (Stratton/Wolfram-without-the-extra-sign). I added a sentence explaining the global chi_n sign cancels in the fractional-linear recursion and cross-products, so both conventions give identical observables. The recursion/cross-product equations are written abstractly and hold under either convention.
- The tau_n Legendre identity is convention-dependent: under scipy lpmv's Condon-Shortley phase, code's 'tay' = -tau_n. I made this explicit (Eq. tau-id gives dP_n^1/dtheta = (1/2)[n(n+1)P_n^0 - P_n^2], and tay computes the negative). Observables are unaffected because pi_n and tau_n carry the same (-1)^1 factor and the pattern is taken in modulus. A reader using a non-CS-phase Legendre routine must drop this relative sign.
- Eq. (psi-deriv) now states the recurrence for chi_n' holds 'by replacing j->y'; strictly chi_n=-x y_n so chi_n' inherits the same minus sign, but since chi enters only via the sign-invariant ratios this is harmless. The code's _psi_chi uses chi=+x y_n consistently, so its chi_p formula is the literal y-analogue of Eq. (psi-deriv).

### Конфокальный сфероид
**Исправлено:**
- Cross-section formulas: corrected coefficient 4pi->2pi, added missing minus sign in Cext, replaced ||T a||^2 with orientation-averaged Tr(T^dag T)  
  _источник:_ Mishchenko, Travis, Lacis 2002: <Cext>=-(2pi/k^2)Re Tr(T), <Csca>=(2pi/k^2)Tr(T^dag T); confirmed via TERMS user guide (arXiv:2111.12556) and SMARTIES (arXiv:1511.00798)
**Открытые оговорки:**
- Cross-section corrected to standard Mishchenko T-matrix convention: extinction now carries the minus sign and the prefactor is 2pi/k0^2 (was erroneously 4pi/k0^2 with no minus); scattering written as the orientation-averaged Frobenius-norm trace Tr(T^dag T) rather than ||T a||^2. A fixed-orientation variant (with 4pi/k0^2 and the incident-amplitude bilinear form) is given as a parenthetical for completeness. The original section conflated the two normalizations.
- Full text of Asano & Yamamoto 1975 and Farafonov & Voshchinnikov 2012 is paywalled; the coupling structure (m-decoupling, l-coupling, TE/TM split at m=0, confocal requirement, EBCM/spheroidal basis) is verified from the published abstracts and from DLMF/Stratton for the special-function machinery, which is sufficient to support every load-bearing claim. The exact per-element form of the boundary block Q^(m) and of alpha/beta projection coefficients is schematic in the section (stated as 'composed of R^(1,3), their derivatives, and M,N combinations') and is presented as such, not as a verbatim quoted formula, so no literal-formula mismatch arises.
- eq:overlap normalization (denominator = self-overlap of S_ml'(c_{j+1})) is a reasonable projection convention but the precise normalization used by Farafonov differs in detail (they project fields/potentials, not bare angular functions); the qualitative statement (non-orthogonal bases at different c_j induce an l-coupling overlap matrix that reduces to delta_ll' when c_j=c_{j+1}) is correct and is what the text asserts.

### Триаксиальный эллипсоид
**Исправлено:**
- Lame harmonics count — section originally said '-n<=m<=n giving 2n+1 functions'; corrected to standard Lame index range 1<=m<=2n+1  
  _источник:_ Standard Lame/ellipsoidal-harmonic indexing (Stratton 1941; DLMF Ch. 29/30): for degree n there are 2n+1 ellipsoidal harmonics indexed m=1..2n+1, not -n..n. Count (2n+1) was right; index range was the spherical-harmonic convention and was corrected.
- eq:two-layer-eff (ORIGINAL effective-permittivity formula) — REPLACED  
  _источник:_ The original eps^eff formula did NOT reproduce the established confocal coated-ellipsoid polarizability (Sihvola 1990; Kostinski & Mongkolsittisilp 2013, JQSRT 131, Eq. 7) under any role assignment of L^(1)/L^(2) — numerically off by 16-38% — and would fail the coated-sphere (L=1/3) limit. Replaced with the directly-verified two-layer polarizability eq:two-layer-alpha, which symbolically reduces EXACTLY to the Bohren-Huffman coated-sphere polarizability (ratio=1) at L^(1)=L^(2)=1/3 and to eq:alpha-homog as f->0.
- eq:two-layer-alpha (NEW) — confocal core-shell ellipsoid polarizability  
  _источник:_ Kostinski & Mongkolsittisilp 2013 (JQSRT 131, 194), Eq. 7, generalized to host eps_m (their Eq. 7 is the eps_m=1 case); Sihvola & Lindell 1990. Verified symbolically to equal Bohren-Huffman coated-sphere (BH Eq. 5.36) when all L=1/3.
**Открытые оговорки:**
- The transfer-matrix Eqs. (eq:transfer-product, eq:transfer-matrix) and the Wronskian normalization are presented at the correct STRUCTURAL level (Sihvola & Lindell 1990 use a transmission-line / 2x2 transfer formalism), but I could not retrieve the verbatim matrix entries from Sihvola & Lindell 1990 (paywalled; arXiv not available). The element-by-element form in eq:transfer-matrix is plausible and dimensionally/sign consistent but is reconstructed rather than copied from the source. The scalar polarizability results eq:two-layer-alpha and eq:alpha-eff ARE rigorously verified (reduce exactly to the coated sphere), so the physics is correct; only the explicit 2x2 entries should be treated as illustrative unless cross-checked against the original Sihvola-Lindell paper.
- DLMF 19.33 states the depolarization/demagnetizing sum as L_a+L_b+L_c=4*pi in the Gaussian (cgs) convention; the section uses the SI/Bohren-Huffman normalization where the sum is 1 (factor 4*pi absorbed). This is a convention difference, not an error — the section is internally consistent with the abc/2 prefactor in eq:depol-integral, which is the SI form. No change needed, but readers comparing to DLMF should note the 4*pi convention.
- The Lame wave equation eq:lame-wave is given in the standard algebraic Jacobi form; the exact value/normalization of the separation constants alpha_0, alpha_1 is left implicit (as is conventional). This is acceptable for a theory overview but is not a fully closed specification.

### Конечный цилиндр (полу-аналитика)
**Исправлено:**
- Far-field amplitude phase factors: M-family (-i)^{n+1}, N-family (-i)^n times e^{ikr}/(kr)  
  _источник:_ Bohren & Huffman (bohren1983); Mishchenko 2002 (mishchenko2002). Draft used a single (-i)^n for BOTH families. Corrected: M-term carries (-i)^{n+1} (from h_n^{(1)}->(-i)^{n+1}e^{ikr}/(kr)), N-term carries (-i)^n (from derivative combination). Standard VSWF far-field asymptotics.
**Открытые оговорки:**
- SCOPE MISMATCH: The task prompt asks to verify the CONE (sphero-conal separation, non-integer Legendre degree nu characteristic equation, Klinkenbusch semi-infinite-cone dyadic Green's function). The file /tmp/sec_cylinder.tex contains ONLY the finite-cylinder section -- there is NO cone content in it. The cone is evidently a separate section in another file. I verified the cylinder section fully; the cone claims could NOT be checked against this file because they are absent here. If a cone section exists, point me at its file (e.g. /tmp/sec_cone.tex) for a separate pass.
- Cone background nonetheless verified against sources for when that section is reviewed: Klinkenbusch 2007 (klinkenbusch2007, Radio Sci. 2007RS003649) treats the semi-infinite circular AND elliptic cone via spherical-multipole analysis in SPHERO-CONAL coordinates, with angular eigenfunctions being periodic and NON-PERIODIC LAME functions (not, in general, associated Legendre P_nu). The circular cone is the special case where Lame functions reduce to non-integer-degree associated Legendre functions P_nu^m; the characteristic equation fixing nu comes from P_nu^m(cos theta_0)=0 (Dirichlet/E-type) or its theta-derivative =0 (Neumann/H-type) on the cone surface theta=theta_0, with the PEC condition combining both. The draft prompt's phrasing 'NON-INTEGER-degree Legendre P_nu for the cone' is correct ONLY for the circular cone; for the general (elliptic) cone the precise framework is Lame functions in sphero-conal coordinates. This nuance should be stated in the cone section.
- PAYWALL/ACCESS: Klinkenbusch 2007 (Wiley/AGU 10.1029/2007RS003649) returned HTTP 403; the Kiel CEM publications page refused connection. Mishchenko 2002 book PDF and several arXiv PDFs returned as undecodable binary to the fetch tool. The cone characteristic equation and the exact dyadic Green's function form (free-space and cone, in sphero-conal coords, in terms of spherical Bessel/Hankel of non-integer degree) were therefore confirmed only at the level of secondary descriptions/abstracts and standard-theory cross-checks, not from the verbatim primary equations.
- Convention note (not an error): the section uses the normalized far-field convention E_sca -> F e^{ikr}/(kr) with dimensionless F and the matching 1/k^2 cross-section prefactors, rather than Bohren & Huffman's dimensional E_sca -> F e^{ikr}/r. Both are self-consistent; ensure the rest of the manuscript (the referenced tmatrix-eq:* and spheroid-eq:cross) uses the SAME normalization, otherwise prefactors will clash across sections.
- Unverified-by-primary-source but standard: the claim that the multilayer radial transfer for k_z != 0 is a 4x4 matrix coupling (m,n) and splits into two 2x2 channels at k_z=0 is physically correct and consistent with the oblique-incidence coupled-wave structure, but no single citation was located that states the 4x4/2x2 split verbatim; it follows from the M/N coupling already cited (Stratton/Bowman).

### Конечный конус (полу-аналитика)
**Исправлено:**
- Cone-tip field behavior: potential ~ r^nu, field components ~ r^{nu-1}; singular field (nu_min<1) only for wide/re-entrant sharp conducting spike, NOT thin needle (CORRECTED from draft which said field~r^nu and 'sharp cone nu<1')  
  _источник:_ Idemen, 'Electric singularity near tip of a sharp cone' (IEEE) / 'Tip singularity at apex of material cone' (Wave Motion): E ~ R^{nu-1}; transcendental eqn in Legendre functions of complex order; thin cone -> regular (consistent with own theta_0->0 limit nu_min->infty)
- Cross sections: C_ext=-(2pi/k0^2)Re[d^dag T d], C_sca=(2pi/k0^2)||Td||^2 (CORRECTED from draft's 4pi/k0^2)  
  _источник:_ Mishchenko 2002 / TERMS guide (arxiv 2111.12556): orientation-avg <C_ext>=-(2pi/k^2)Re sum T_diag, <C_sca>=(2pi/k^2) sum|T|^2; fixed-orientation cross sections proportional to same right-hand sides with 2pi/k^2 prefactor in power-normalized VSWF convention
**Открытые оговорки:**
- Klinkenbusch 2007 (Radio Sci. 10.1029/2007RS003649) is paywalled (HTTP 403); verification of its dyadic Green's function and boundary-condition details relied on the journal abstract, the open URSI procGA08 B02p7 elliptic-cone paper, and the arxiv 2512.18498 spherical-cavity-with-cone paper rather than the primary text. The schematic Green's function eq. (cone-eq:greendyad) and the M/N grouping match standard form but were not checked line-by-line against Klinkenbusch's exact normalization N_nu.
- Cross-section prefactor is convention-dependent. Corrected to 2pi/k0^2 to match the Mishchenko power-normalized VSWF convention (cited mishchenko2002) and the orientation-averaged forms in the TERMS guide. If the document's plan section (plan-eq:Tdef) uses a different VSWF normalization (e.g. unnormalized Stratton VSWFs), the prefactor and even the form ||Td||^2 would change; the editor should confirm consistency with the global normalization fixed in the sphere/plan sections.
- Cone-tip exponent corrected to field ~ r^{nu_min-1} (Idemen) and singular case (nu_min<1) attributed to a wide/re-entrant conducting spike rather than a thin needle, removing the internal contradiction with the theta_0->0 limit (nu_min->infty, regular). The precise mapping of which physical PEC configuration yields nu_min<1 depends on the document's interior/exterior convention (line 23 has it ambiguous: exterior 0<=theta<=theta_0); the editor should make the scatterer-fills-which-region convention explicit and consistent across cone-eq:eig and the theta_0->pi limit bullet.
- TM<->N(electric) and TE<->M(magnetic) Debye-potential assignment with Dirichlet on the electric potential and Neumann on the magnetic potential is standard and consistent with search results, but the exact pairing (which potential's BC is Dirichlet vs Neumann) could not be confirmed verbatim from a primary open source; it agrees with the general statement 'Dirichlet on one potential, Neumann on the conjugate' but is asserted, not quoted.
- Aden-Kerker recursion eq. (cone-eq:adenkerker) cross-references the document's own sphere section (sphere-sec:sphere) and could not be checked against the allowed external citation keys; its non-integer-nu generalization is plausible but unverified against primary literature.
- Meixner exponent tau is stated to 'depend on the dihedral angle of the rim' without an explicit value; this is correct in spirit (edge exponent depends on wedge angle) but no closed-form was verified for the specific cone-base dihedral; left as qualitative.

### Суперпозиция T-матриц (GMM)
**Исправлено:**
- Cluster T-matrix via regular translations Rg-A,Rg-B about common origin (eq tm-cluster)  
  _источник:_ Mishchenko et al. 2002, Eqs. (5.232)-(5.236): cluster T uses RgA/RgB (regular, j_l) operators; corrected argument labels to r_{jO}, r_{Oi} consistent with (5.236)
- Scattering cross section C_sca = (1/k^2) int |F|^2 dOmega = (1/k^2) sum |c|^2 (eq tm-Csca)  
  _источник:_ Mishchenko et al. 2002, Eq. (5.18b): C_sca = (1/(k^2 |E|^2)) sum (|p|^2+|q|^2); replaced bare (...) placeholder with explicit overlap-integral form reducing to sum of squares
- Extinction cross section sign and form (eq tm-Cext)  
  _источник:_ Mishchenko et al. 2002, Eq. (5.18a): C_ext = -(1/(k^2|E|^2)) Re sum[a_mn p_mn* + b_mn q_mn*] — leading MINUS with e^{-iwt}; removed the self-contradictory '-(4pi/k^2)Re[...] = +(4pi/k^2)Re sum(...)' double-sign placeholder
**Открытые оговорки:**
- Cruzan A/B coefficient formulas (eq tm-Acoef/tm-Bcoef): the exact prefactor and the precise definition/normalization of the Gaunt coefficients a(m,n,-mu,nu,l) and b(...) could not be byte-verified against Cruzan 1962 because the publisher PDF (AMS/DTIC/scispace) returned HTTP 403/binary. The structure (i^{nu-n}(2n+1)/(2n(n+1)), factorial ratio, [nu(nu+1)+n(n+1)-l(l+1)] weight for A, sum over l on Gaunt selection rules, Y_l^{m-mu}) is the standard Cruzan/Xu form and is internally consistent; the text already defers exact normalization to the primary sources (Cruzan 1962, Stein 1961, Xu 1995). A reader implementing should take a,b from Xu 1995 recurrences as the section recommends.
- The 4pi vs 1/k^2 prefactor equality in eq tm-Cext (coefficient form -(1/k^2)Re sum[d c*] vs forward-amplitude form (4pi/k^2)Re[e*.F]) depends on the precise normalization of the angular functions X_mn and of d_mn relative to F; both forms are individually correct in Mishchenko's convention (Eqs. 5.18a and 5.97), but the literal numerical equality between the two written expressions holds only with Mishchenko's specific C_mn/B_mn normalization. This is flagged as a convention-coupling, not an error.
- The angular overlap integrals I^{MM}, I^{NN} in eq tm-Csca are written schematically (with '+...' for the inter-type/cross terms); the full closed form (Xu 1995 'far field' paper, Appl. Opt. 36, 9496) was not transcribed. For single-body L=1 they correctly collapse to sum|c|^2 (Mishchenko Eq. 5.18b).

### Разложение и верификация
**Исправлено:**
- Clausius-Mossotti reduction at L_j=1/3  
  _источник:_ Setting L_j=1/3 in eq:depol gives alpha=3V(eps-eps_m)/(eps+2eps_m)=4 pi a^3 (eps-eps_m)/(eps+2eps_m). The original LaTeX had a SPURIOUS extra factor eps_m: '3V eps_m (eps-eps_m)/(eps+2eps_m)'. Removed it. Confirmed: standard sphere polarizability alpha=4 pi a^3 (eps-eps_m)/(eps+2eps_m), Bohren & Huffman / Clausius-Mossotti (farside.ph.utexas.edu, kirkmcd.princeton.edu).
**Открытые оговорки:**
- Far-field/cross-section: the section states the optical-theorem cross-section formula but does not give the explicit Q_sca = (2/x^2) sum (2n+1)(|a_n|^2+|b_n|^2) series (it lives in analytic_mie.py). This is acceptable since the section's scope is decomposition+verification, not the per-primitive far-field derivation, which presumably appears in the sphere section. Not an error.
- The Riccati relation xi_n = psi_n - i chi_n in the verification preamble is consistent with the e^{-iwt}/h_n^(1) convention (h_n^(1)=j_n + i y_n, chi_n = -x y_n => xi_n = x j_n - i(-x y_n)? note xi_n = x h_n^(1) = x j_n + i x y_n = psi_n - i chi_n with chi_n := -x y_n). The notation is the standard Bohren-Huffman one and is internally consistent; no change needed.
- SMARTIES is cited via mishchenko2002 as a stand-in reference since no dedicated SMARTIES bib key exists in the allowed list; the citation supports T-matrix spheroid methodology generally but is not the SMARTIES paper itself. This is a minor citation-precision compromise forced by the fixed bib key list, not a factual error.
