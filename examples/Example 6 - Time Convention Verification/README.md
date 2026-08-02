# Example 6 — Time Convention Verification (e^{+jωt})

Verifies that GreenTensor's scattering coefficients `Mn`, `Nn`
(`green_tensor/calc.py`, class `RCSCalculator`) follow the canonical
engineering time convention **e^{+jωt}** adopted throughout the theory:

* outgoing radial dependence — Riccati–Hankel functions of the **second**
  kind, ζₙ(z) = z·hₙ⁽²⁾(z) = ψₙ(z) + j·χₙ(z);
* lossy media — ε = ε′ − jε″.

## How the convention is implemented in the code

`scipy.special.hankel1` is native to the opposite (physics, e^{−iωt})
convention. `calc.py` therefore evaluates the coefficient ratio with
h⁽¹⁾ and then applies an explicit conjugation:

```python
# green_tensor/calc.py, calculate_scattering_coefficients()
self.Mn[i] = (Z*mJ - mJpr) / (Z*mH - mHpr)   # computed with hankel1
self.Mn[i] = self.Mn[i].real - self.Mn[i].imag * 1j   # -> e^{+jωt} canon
```

For real arguments `conj(h⁽¹⁾) = h⁽²⁾` while ψ, χ and the propagated
impedances Z, Y are real, so the conjugation is mathematically identical
to evaluating the ratio with h⁽²⁾ directly. This example **proves that
identity numerically** against an independent solver written directly in
the e^{+jωt} convention (ζ = ψ + jχ, no conjugation anywhere).

## Checks

| # | Check | Result |
|---|-------|--------|
| V1 | Single-layer sphere (ε = 2.5, k₀a = 2π): `Mn`,`Nn` vs direct-h⁽²⁾ aₙ,bₙ | ≈ 7·10⁻¹⁵ |
| V2 | 8-layer Fuchs-optimized Luneburg lens (k₀a = 4π): same identity | ≈ 3·10⁻¹⁴ |
| V3 | Outgoing wave: d(arg ζₙ(k₀r))/dr → −k₀ (phase fronts move outward with e^{+jωt}) | −0.9994 |
| V4 | Energy (lossless): Re Mₙ = \|Mₙ\|² (unitarity circle) and Q_ext = Q_sca (optical theorem) | ≈ 5·10⁻¹⁵ |
| V5 | Observable invariance: backscatter σ_r(k₀a) sweep, GreenTensor vs direct h⁽²⁾ | ≈ 6·10⁻¹⁴ |

## Run

```bash
python3 time_convention_verification.py
```

Produces `time_convention_verification.png/.pdf`:
left — Mₙ, Nₙ of the 8-layer lens on the unitarity circle, both solvers
overlaid; right — backscatter σ_r/(πa²) of the same lens over
k₀a ∈ [1, 4π], line (GreenTensor) vs markers (direct h⁽²⁾).

## Reference

Dissertation, ch. "Tensor Green's function", eq. (Mₙ, Nₙ from Z, Y) —
the h⁽²⁾ form of the scattering coefficients; `GreenTensor_Theory.pdf`.
