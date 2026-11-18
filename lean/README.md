# Lean formalisation

Lean 4.35.0 with Mathlib v4.35.0. The development uses no `sorry`, no custom `axiom` and no
`native_decide`. For every listed theorem, `scripts/print_axioms.lean` reports only `propext`,
`Classical.choice` and `Quot.sound`.

There are two ways to build:
- **Standalone:** `lake exe cache get && lake build`.
- **Script:** `scripts/build_with_lean_path.sh` builds with plain `lean` via `LEAN_PATH` and prints the
  axioms; `scripts/verify.py` uses it.
  - By default it uses this directory's own Lake environment, `.lake/packages`, so run
    `lake exe cache get` first.
  - To reuse an existing Mathlib v4.35.0 build instead, set `MATHLIB_PACKAGES` to its `.lake/packages`
    directory. This is an explicit opt-in.

## What is formalised

| File | Theorem(s) | Content |
|---|---|---|
| `AreaCorollary.lean` | `amgm3`, `amgm3_eq` | Three-variable AM–GM, 27xyz ≤ (x+y+z)³, with its equality case |
| | `triangle_isoperimetric_heron`, `triangle_isoperimetric_heron_eq` | 27·(a+b+c)(−a+b+c)(a−b+c)(a+b−c) ≤ (a+b+c)⁴ for a nondegenerate triangle, with equality only if a = b = c. With Heron's formula 16·Area² = (product), this gives Area/L² ≤ √3/36 |
| | `area_corollary`, `area_corollary_strict` | If σL ≤ S (resp. < S) and A ≤ kL², then σ√A ≤ S√k (resp. <). This is Corollary 1.2 |
| `Reformulation.lean` | `G_spec` | (AᵀA)·G = I for the stated entries of G = A⁻¹A⁻ᵀ, with A = [[1, ρa], [0, 1+ρb]] |
| | `X_eq_rho_Y` | Entrywise identity X = ρY of Lemma 4.2, from Σ_e μᵉ = 0 and Q₁₁ + Q₂₂ = σ Σ_e gᵉ. In Lean, the paper's μᵉ are `m0x`…`m2y`, ωₑ are `c0`, `c1`, `c2`, λₑ are `d1`, `d2`, and 𝒩_A, 𝒟_A are `NA`, `DA` |
| | `stretch_div`, `stretch_sq1`, `stretch_sq2` | (s−1)/ρ = (κ+ρ)/(1+s) from s² = 1 + ρκ + ρ², and the two values of κ |
| `CosBound.lean` | `int_x_sin_sq`, `int_cos_sq`, `int_cos` | The closed-form integrals of cos(πx). The first gives 𝓘 in Lemma 5.1; the other two are used only in the consistency checks: ∫₀ᵖ x sin²(πx), ∫₀ᵖ cos²(πx), ∫₀ᵖ cos(πx) |

## What is assumed, not formalised

- Heron's formula.
- The reduction of the triangle integrals in Lemma 5.1 to one-dimensional integrals (Fubini over
  vertical sections, and the arc-length parametrisation of the sides).
- The whole spectral theory: the variational principle, the DtN map, the inclusion lemma, Davis–Kahan,
  Schur's lemma, the CR projection bound, the trace and Payne–Weinberger inequalities.
- All ball-arithmetic computations, which are done in Python/Arb.
