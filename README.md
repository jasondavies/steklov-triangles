# The equilateral triangle maximises the first Steklov eigenvalue among triangles

This repository contains a paper and its verification material. It answers the case n = 3 of the Steklov
Pólya–Szegő question (AIM problem list *Steklov eigenproblems*, Problem 4.6, É. Oudet): among triangles,
the equilateral triangle E uniquely maximises the first nonzero Steklov eigenvalue σ₁, both at fixed
perimeter and at fixed area. The proof is computer-assisted.

## Results

Write Σ(T) = σ₁(T)·L(T), where L is the perimeter.

- **Theorem 1.1.** For every triangle T, Σ(T) ≤ Σ(E), with equality only for equilateral T, and
  3.8724615957722 ≤ Σ(E) ≤ 3.8724615957734.
- **Corollary 1.2 (fixed area).** σ₁(T)·√|T| ≤ σ₁(E)·√|E| ∈ [0.8494076787253, 0.8494076787257], with
  equality only for equilateral T.
- **Quantitative stability.** Normalise a longest side to a unit base and let ρ be the distance of
  the apex from the equilateral apex. The deficit Σ(E) − Σ(T) is at least ρ/2, hence at least
  (ℓmax − ℓmin)/(2ℓmax). A deficit of 0.01 or less forces ℓmin/ℓmax ≥ 0.98.
- **Equilateral spectrum.** For side length one, the cubic eigenfunction spans the eigenspace of
  σ₃(E) = 2√3, and σ₄(E) > 3.5617858371486.
- **Splitting of the first pair.** Displacing the apex by ρ over a fixed unit base gives first-order
  changes −cρ and +cρ in σ₁L and σ₂L, with uniform O(ρ²) remainders and
  3.5332113353 < c < 3.5332113355.
- **Thin triangles.** For each fixed apex abscissa 0 < p < 1 and height q tending to zero,
  Σ(T) = 2qθ(p) + o(q), by the thin-domain lemma of Henrot and Michetti.

Cheng–Gui–Hu–Li–Yao (arXiv:2603.25116v2, 2026) record the conjecture as open for every N ≥ 3. We found no
earlier proof of the triangle case.

## Read

- [Manuscript](paper/steklov-triangles.pdf) and its [LaTeX source](paper/steklov-triangles.tex).
- [code/README.md](code/README.md): the scripts.
- [lean/README.md](lean/README.md): the Lean development.

## Reproduce

You need:
- Python 3.12 with the pins in `requirements.txt`: python-flint 0.9.0, numpy 2.5.5, scipy 1.18.1,
  mpmath 1.4.1;
- `latexmk`/`pdflatex`;
- Lean 4.35.0 with Mathlib v4.35.0, obtained by `cd lean && lake exe cache get`.

```sh
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
(cd lean && lake exe cache get)                     # Mathlib for the Lean step
.venv/bin/python scripts/verify.py                  # everything (a few minutes on one core)
.venv/bin/python scripts/verify.py --only certified|checks|numerics|paper|lean
.venv/bin/python scripts/verify.py --skip-numerics  # skip the float-only numerics
MATHLIB_PACKAGES=/path/to/.lake/packages .venv/bin/python scripts/verify.py   # reuse an existing Mathlib build
```

`scripts/verify.py` works as follows:
- Each step writes a JSON file to `verification/`.
- Each gate compares that output with the paper's claims in [claims.json](claims.json), in exact
  rational arithmetic.
- The paper's numbers come from the same file, through the generated macros `paper/claims.tex`
  (`scripts/make_claims_tex.py`).
- Logs go to `verification/<step>.log`, and a summary with input and PDF hashes to
  `verification/summary.json`.
- It exits nonzero on any failed gate, LaTeX warning or overfull box. It uses one core (BLAS threads = 1).

**Certificate data** ([data/](data/)) is committed:
- `equilateral_pair.json`: the binary64 coefficients c_k of the equilateral eigenpair approximation;
- `cr_congruence.npz`: the congruence matrix of the Crouzeix–Raviart inertia certificate;
- `table_estimates.json`: the high-precision table data and its trial functions.

The certified steps read only this data, so their output does not depend on the BLAS library or the CPU.
Two optional scripts regenerate the data; neither is run by `verify.py`:
- `code/make_certificate.py` for the first two files;
- `code/hp_table.py` for the table data, which is slow.

## What is proved, certified, formalised, or only numerical

**Proved by hand in the paper:**
- the reduction to the parameter domain;
- the inclusion and eigenvector-enclosure lemmas for the Dirichlet-to-Neumann map;
- the Crouzeix–Raviart trace constant and projection bound;
- the cubic boundary condition, spectral counting and D₃ multiplicity argument;
- the affine transplantation and the factorisation X = ρY;
- the description of the far region and the cos(πx) bound;
- the uniform first-order splitting via the cluster projection and cyclic relabelling;
- the stability and fixed-area inequalities.

The cited inputs are:
- the Payne–Weinberger inequality;
- standard facts about the Dirichlet-to-Neumann map;
- Gershgorin's theorem and Sylvester's law of inertia;
- Henrot and Michetti’s Lemma 3.2 for the fixed-profile thin limit.

**Certified in Arb ball arithmetic** (python-flint, `code/`):
- **Σ(E).** A harmonic-polynomial pair with relative residual η ≤ 1.96·10⁻¹³ gives the enclosure of Σ(E)
  (`equilateral.py`).
- **σ₄(E) > 3.5617858371486.** Crouzeix–Raviart elements with 900 dofs and an inertia certificate
  (`cr_lower.py`).
- **Local region, ρ ≤ ρ₁.** 907 boxes in (ρ,ψ) (`local.py`). They are re-checked, and their exact tiling
  of the parameter rectangle is verified, by `coverage_check.py`. Both regions use the nearest binary64 value to 0.6 for ρ₁.
- **Far region, ρ ≥ ρ₁.** 16 intervals in p, with upper bound 3.114382798607 (`far_region.py`).
- **Table of Section 7.** A rigorous upper bound for Σ of each listed polygon, from a fixed trial
  function on the exact polygon (`table_bounds.py`). The estimated values themselves are numerical.

**Formalised in Lean 4/Mathlib** (no `sorry`; axioms only `propext`, `Classical.choice`, `Quot.sound`):
- the triangle isoperimetric inequality in Heron form, with its equality case;
- the area corollary;
- the algebra of X = ρY;
- the closed-form integrals of cos(πx).

**Numerical only.** The Ritz maps, estimated table errors, and sampled values of the one-dimensional thin-profile eigenvalue.

**Main residual risks.** The by-hand derivations of the Crouzeix–Raviart constant and of the closed forms,
their transcription into code, and the correctness of python-flint/Arb. No independent human review has
been recorded.

## Authorship and AI disclosure

Author: Jason Davies. The mathematics, code, Lean development and manuscript were produced by AI systems
(Anthropic’s Claude and OpenAI Codex) under the author's direction. The author has not checked the arguments; see the
declaration in the paper.

## License

Apache License 2.0; see [LICENSE](LICENSE). Third-party dependencies keep their own licences.
