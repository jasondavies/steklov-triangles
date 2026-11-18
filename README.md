# The equilateral triangle maximises the first Steklov eigenvalue among triangles

This repository contains a paper and its verification material. It answers the case n = 3 of the Steklov
Pólya–Szegő question (AIM problem list *Steklov eigenproblems*, Problem 4.6, É. Oudet): among triangles,
the equilateral triangle E uniquely maximises the first nonzero Steklov eigenvalue σ₁, both at fixed
perimeter and at fixed area. The proof is computer-assisted.

## Results

Write Σ(T) = σ₁(T)·L(T), where L is the perimeter.

- **Theorem 1.1.** For every triangle T, Σ(T) ≤ Σ(E), with equality only for equilateral T. Also
  σ₁(E) = σ₂(E) < σ₃(E), and 3.8724615957722 ≤ Σ(E) ≤ 3.8724615957734.
- **Corollary 1.2 (fixed area).** σ₁(T)·√|T| ≤ σ₁(E)·√|E| ∈ [0.8494076787253, 0.8494076787257], with
  equality only for equilateral T.
- **Numerical observations** (not rigorous, except where marked certified):
  - a table of Σ for seven triangles and the regular 4-, 5- and 6-gons. Each entry gives a certified Arb
    upper bound and a high-precision estimate with an estimated error;
  - near E, Σ ≈ Σ(E) − 3.533ρ, so the maximum is a cone;
  - for thin triangles, Σ ≈ 2q·θ(p);
  - σ₃(E) agrees with 2√3 to 14 digits. 2√3 is an eigenvalue with an explicit eigenfunction, but equality
    σ₃(E) = 2√3 is not proved.

Cheng–Gui–Hu–Li–Yao (arXiv:2603.25116v2, 2026) record the conjecture as open for every N ≥ 3. We found no
earlier proof of the triangle case.

## Read

- [Manuscript](paper/steklov-triangles.pdf) and its [LaTeX source](paper/steklov-triangles.tex), 11 pages.
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
Runs with `OPENBLAS_CORETYPE=Sandybridge` and `=Haswell` produce identical certificates. Two optional
scripts regenerate the data; neither is run by `verify.py`:
- `code/make_certificate.py` for the first two files;
- `code/hp_table.py` for the table data, which is slow.

## What is proved, certified, formalised, or only numerical

**Proved by hand in the paper:**
- the reduction to the parameter domain;
- the inclusion and eigenvector-enclosure lemmas for the Dirichlet-to-Neumann map;
- the Crouzeix–Raviart trace constant and projection bound;
- the D₃/Schur multiplicity argument;
- the affine transplantation and the factorisation X = ρY;
- the description of the far region and the cos(πx) bound;
- the fixed-area corollary.

The cited inputs are:
- the Payne–Weinberger inequality;
- standard facts about the Dirichlet-to-Neumann map;
- Gershgorin's theorem and Sylvester's law of inertia.

**Certified in Arb ball arithmetic** (python-flint, `code/`):
- **Σ(E).** A harmonic-polynomial pair with relative residual η ≤ 1.96·10⁻¹³ gives the enclosure of Σ(E)
  (`equilateral.py`).
- **σ₃(E) ≥ 2.2003.** Crouzeix–Raviart elements with 408 dofs and an inertia certificate
  (`cr_lower.py`).
- **Local region, ρ ≤ 0.6.** 460 boxes in (ρ,ψ) (`local.py`). They are re-checked, and their exact tiling
  of the disc is verified, by `coverage_check.py`.
- **Far region, ρ ≥ 0.6.** 16 intervals in p, with largest bound 3.114 < 3.8725 (`far_region.py`).
- **Table of Section 7.** A rigorous upper bound for Σ of each listed polygon, from a fixed trial
  function on the exact polygon (`table_bounds.py`). The estimated values themselves are numerical.

**Formalised in Lean 4/Mathlib** (no `sorry`; axioms only `propext`, `Classical.choice`, `Quot.sound`):
- the triangle isoperimetric inequality in Heron form, with its equality case;
- the area corollary;
- the algebra of X = ρY;
- the closed-form integrals of cos(πx).

**Numerical only.** The Ritz maps and special values, the cone slope, the thin-triangle asymptotics, and
σ₃(E) ≈ 2√3. These checks also confirm that every certified bound lies above the converged Ritz value.

**Main residual risks.** The by-hand derivations of the Crouzeix–Raviart constant and of the closed forms,
their transcription into code, and the correctness of python-flint/Arb. No independent human review has
been recorded.

## Authorship and AI disclosure

Author: Jason Davies. The mathematics, code, Lean development and manuscript were produced by AI systems
(Anthropic's Claude) under the author's direction. The author has not checked the arguments; see the
declaration in the paper.

## License

Apache License 2.0; see [LICENSE](LICENSE). Third-party dependencies keep their own licences.
