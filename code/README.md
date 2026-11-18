# Computational material

Run everything with `.venv/bin/python scripts/verify.py` from the repository root. Run individual
scripts from `code/`. Every step writes a machine-readable JSON file to `../verification/`. In those
files, certified endpoints are exact `man*2^exp` numbers; `verify.py` compares them with
`../claims.json`.

The certified steps read only the committed certificate data in `../data/`:
- `equilateral_pair.json`: the coefficients c_k;
- `cr_congruence.npz`: the congruence matrix O.

They therefore do not depend on BLAS or floating-point eigensolvers. `make_certificate.py` regenerates
this data. It is optional and is not run by `verify.py`.

**Certified** (Arb ball arithmetic via python-flint; together these make up the proof):

| Script | What it certifies | Paper |
|---|---|---|
| `rig.py` | Library: exact polynomial line integrals along edges; Green-formula energy moments; edge Gram matrices and means; residual Gram matrices | §3, §4 |
| `equilateral.py` | Reads the coefficients c_k and builds the exactly D₃-equivariant harmonic pair (degree 20). Certifies the relative residual η ≤ 1.9559e-13 and the enclosures of Σ(E) and σ₁√\|E\|. Writes `equilateral.json` | §3.1, Prop. 3.5 |
| `cr_lower.py` | Crouzeix–Raviart mesh (n_h = 16) and inertia of S_h − 3B_h, with exactly 3 negative eigenvalues, using the committed congruence matrix. Gives σ₃(E) ≥ 2.2003. Writes `cr_lower.json` | Lemmas 3.1–3.2, Prop. 3.3 |
| `local.py [rho1] [nψ]` | Moment enclosures of the exact eigenspace and the box covering of [0, ρ₁] × [0, 6.3]; writes `local_leaves.json` | Lemmas 4.1–4.3, Prop. 4.4 |
| `coverage_check.py` | The local boxes tile [0, ρ₁] × [0, 6.3] (exact area sum, disjoint interiors) and are re-certified; ρ₁ is the same double as in `far_region.py` | Prop. 4.4 |
| `far_region.py` | One-dimensional bound for the cos(πx) trial function on 16 p-intervals covering the far region | Lemma 5.1, Prop. 5.2 |

**Consistency checks** (floating point; not part of the proof):

| Script | Check |
|---|---|
| `test_moments.py` | Green-formula moments vs brute-force quadrature |
| `test_independent.py` | Full cos(πx) Rayleigh quotient and transplant bound vs direct quadrature on T |
| `check_local_vs_ritz.py` | Every bound dominates the floating-point Ritz value |
| `table_bounds.py` | Certified Arb upper bounds for the table of §7: the Rayleigh quotient of the committed trial functions on the exact polygons |

**Numerics** (floating point, not rigorous):

| Script | Output |
|---|---|
| `hritz.py` | Float64 harmonic-polynomial Ritz solver. It drops ill-conditioned directions, so at high degree its values can be slightly too high |
| `hp_table.py` | High-precision Ritz values (Arb matrices, mpmath eigensolver; degrees 30/40/50, two precisions) and trial functions for the table. Writes `../data/table_estimates.json`. Optional and slow; not run by `verify.py` |
| `make_certificate.py` | Optional: regenerates `../data/equilateral_pair.json` and `../data/cr_congruence.npz` |
| `map_sigma.py` | Σ over a grid on the parameter domain (`map_sigma.json`) |
| `numerics_extra.py` | Thin-triangle asymptotics; the explicit eigenvalue 2√3 of E; the CR discrete values |

Some identifiers differ from the paper's notation:
- the paper's σ̃ is `sig`;
- the edge means μᵉ are `mm` (from `m` in `rig.moments`);
- λ₁, λ₂ are `d1`, `d2`;
- the ωₑ are `ce`;
- 𝒟_A is `D`;
- the function 𝓘 of Lemma 5.1 is `j`.

Dependencies are pinned in `../requirements.txt`. The scripts read no external data.
