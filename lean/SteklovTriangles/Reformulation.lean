import Mathlib

/-!
# The singularity-free reformulation `X = ρ Y` (paper, Lemma 4.2)

All quantities are real numbers; one fixed matrix entry `(x,y)` of the 2×2 matrices is treated
at a time (the identity is entrywise).

* `A = [[1, ρ a], [0, d]]`, `d = 1 + ρ b`;  `G = A⁻¹A⁻ᵀ` has entries
  `G11 = 1 + ρ²a²/d²`, `G12 = G21 = -ρ a/d²`, `G22 = 1/d²`  (`G_spec`: `(AᵀA) G = I`).
* Edge stretches `s0 = 1`, `s1 = 1 + ρ d1`, `s2 = 1 + ρ d2`, perimeter `L = 1 + s1 + s2`.
* Exact eigenspace identities (hypotheses): `m0 + m1 + m2 = 0` for both indices, and the eigen
  relation `Q11 + Q22 = σ (g0 + g1 + g2)` (that is, `N_I = σ D_I`; the mean term of `D_I` vanishes).
* `X = σ (3/L) D_A - N_A` with `D_A = Σ s_e g_e - (Σ s_e m_e^x)(Σ s_e m_e^y)/L`,
  `N_A = d Σ G_ij Q_ij`.
* `Y = σ (Σ c_e g_e - 3 ρ v^x v^y / L²) - Σ H_ij Q_ij` with
  `c0 = -(d1+d2)/L`, `c1 = (2 d1 - d2)/L`, `c2 = (2 d2 - d1)/L`, `v = d1 m1 + d2 m2`,
  `H = [[b + ρ a²/d, -a/d], [-a/d, -b/d]]`.

`X_eq_rho_Y` proves `X = ρ Y`.  `stretch_div` proves `(s-1)/ρ = (κ+ρ)/(1+s)` from
`s² = 1 + ρκ + ρ²`, and `stretch_sq1`, `stretch_sq2` give `κ` for the two moving edges.
-/

namespace SteklovTriangles

theorem G_spec (ρ a b : ℝ) (hd : 1 + ρ * b ≠ 0) :
    let d := 1 + ρ * b
    let G11 := 1 + ρ ^ 2 * a ^ 2 / d ^ 2
    let G12 := -(ρ * a) / d ^ 2
    let G22 := 1 / d ^ 2
    -- AᵀA = [[1, ρ a], [ρ a, ρ²a² + d²]]
    (1 * G11 + (ρ * a) * G12 = 1) ∧ (1 * G12 + (ρ * a) * G22 = 0) ∧
    ((ρ * a) * G11 + (ρ ^ 2 * a ^ 2 + d ^ 2) * G12 = 0) ∧
    ((ρ * a) * G12 + (ρ ^ 2 * a ^ 2 + d ^ 2) * G22 = 1) := by
  intro d G11 G12 G22
  have hd' : d ≠ 0 := hd
  refine ⟨?_, ?_, ?_, ?_⟩ <;> (simp only [G11, G12, G22]; field_simp; ring)

theorem X_eq_rho_Y (σ ρ a b d1 d2 g0 g1 g2 m0x m1x m2x m0y m1y m2y Q11 Q12 Q21 Q22 : ℝ)
    (hd : 1 + ρ * b ≠ 0) (hL : 3 + ρ * d1 + ρ * d2 ≠ 0)
    (hmx : m0x + m1x + m2x = 0) (hmy : m0y + m1y + m2y = 0)
    (heig : Q11 + Q22 = σ * (g0 + g1 + g2)) :
    let d := 1 + ρ * b
    let s1 := 1 + ρ * d1
    let s2 := 1 + ρ * d2
    let L := 1 + s1 + s2
    let G11 := 1 + ρ ^ 2 * a ^ 2 / d ^ 2
    let G12 := -(ρ * a) / d ^ 2
    let G22 := 1 / d ^ 2
    let DA := (g0 + s1 * g1 + s2 * g2) - (m0x + s1 * m1x + s2 * m2x) * (m0y + s1 * m1y + s2 * m2y) / L
    let NA := d * (G11 * Q11 + G12 * Q12 + G12 * Q21 + G22 * Q22)
    let X := σ * (3 / L) * DA - NA
    let c0 := -(d1 + d2) / L
    let c1 := (2 * d1 - d2) / L
    let c2 := (2 * d2 - d1) / L
    let vx := d1 * m1x + d2 * m2x
    let vy := d1 * m1y + d2 * m2y
    let H11 := b + ρ * a ^ 2 / d
    let H12 := -a / d
    let H22 := -b / d
    let Y := σ * (c0 * g0 + c1 * g1 + c2 * g2 - 3 * ρ * vx * vy / L ^ 2)
              - (H11 * Q11 + H12 * Q12 + H12 * Q21 + H22 * Q22)
    X = ρ * Y := by
  intro d s1 s2 L G11 G12 G22 DA NA X c0 c1 c2 vx vy H11 H12 H22 Y
  have hd' : d ≠ 0 := hd
  have hL' : L ≠ 0 := by
    show 1 + (1 + ρ * d1) + (1 + ρ * d2) ≠ 0
    have : 1 + (1 + ρ * d1) + (1 + ρ * d2) = 3 + ρ * d1 + ρ * d2 := by ring
    rw [this]; exact hL
  -- the means: Σ s_e m_e = ρ v  (uses Σ m_e = 0)
  have hsx : m0x + s1 * m1x + s2 * m2x = ρ * vx := by
    simp only [s1, s2, vx]; linear_combination hmx
  have hsy : m0y + s1 * m1y + s2 * m2y = ρ * vy := by
    simp only [s1, s2, vy]; linear_combination hmy
  have hQ : Q22 = σ * (g0 + g1 + g2) - Q11 := by linarith
  simp only [X, DA, NA, Y]
  rw [hsx, hsy]
  simp only [G11, G12, G22, H11, H12, H22, c0, c1, c2]
  rw [hQ]
  simp only [L, s1, s2] at hL' ⊢
  field_simp
  ring

/-- `(s - 1)/ρ = (κ + ρ)/(1 + s)` whenever `s² = 1 + ρκ + ρ²`, `s ≥ 0`, `ρ ≠ 0`. -/
theorem stretch_div (s ρ κ : ℝ) (hs : 0 ≤ s) (hρ : ρ ≠ 0) (h : s ^ 2 = 1 + ρ * κ + ρ ^ 2) :
    (s - 1) / ρ = (κ + ρ) / (1 + s) := by
  have h1 : 1 + s ≠ 0 := by linarith
  rw [div_eq_div_iff hρ h1]
  nlinarith [h]

/-- Edge from `(1,0)` to the apex `(1/2 + ρ cos ψ, √3/2 + ρ sin ψ)`: `κ₁ = -cos ψ + √3 sin ψ`. -/
theorem stretch_sq1 (ρ c s : ℝ) (hcs : c ^ 2 + s ^ 2 = 1) :
    (-(1 : ℝ) / 2 + ρ * c) ^ 2 + (Real.sqrt 3 / 2 + ρ * s) ^ 2
      = 1 + ρ * (-c + Real.sqrt 3 * s) + ρ ^ 2 := by
  have h3 : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  linear_combination (ρ ^ 2) * hcs + (1 / 4 : ℝ) * h3

/-- Edge from the apex to `(0,0)`: `κ₂ = cos ψ + √3 sin ψ`. -/
theorem stretch_sq2 (ρ c s : ℝ) (hcs : c ^ 2 + s ^ 2 = 1) :
    ((1 : ℝ) / 2 + ρ * c) ^ 2 + (Real.sqrt 3 / 2 + ρ * s) ^ 2
      = 1 + ρ * (c + Real.sqrt 3 * s) + ρ ^ 2 := by
  have h3 : Real.sqrt 3 ^ 2 = 3 := Real.sq_sqrt (by norm_num)
  linear_combination (ρ ^ 2) * hcs + (1 / 4 : ℝ) * h3

end SteklovTriangles
