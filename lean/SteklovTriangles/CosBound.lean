import Mathlib

/-!
# Closed-form one-dimensional integrals behind the `cos(πx)` bound (paper, Lemma 5.1)

With `u(x,y) = cos(πx)` on `T(p,q)`:
* `int_x_sin_sq`:  `∫₀ᵖ x sin²(πx) dx = p²/4 - p sin(2πp)/(4π) - (cos(2πp) - 1)/(8π²)`,
  so `𝓘(p) = (1/p) ∫₀ᵖ x sin²(πx) dx` is the closed form used in `code/far_region.py` (`j`).
* `int_cos_sq`:    `∫₀ᵖ cos²(πx) dx = p/2 + sin(2πp)/(4π)`.
* `int_cos`:       `∫₀ᵖ cos(πx) dx = sin(πp)/π`.
The reduction of the triangle integrals to these (Fubini over vertical segments; arc-length
parametrisation of the sides) is elementary and done in the paper, not in Lean.
-/

namespace SteklovTriangles
open Real intervalIntegral

theorem int_x_sin_sq (p : ℝ) :
    ∫ x in (0:ℝ)..p, x * sin (π * x) ^ 2
      = p ^ 2 / 4 - p * sin (2 * π * p) / (4 * π) - (cos (2 * π * p) - 1) / (8 * π ^ 2) := by
  have hπ : π ≠ 0 := pi_ne_zero
  let F : ℝ → ℝ := fun x => x ^ 2 / 4 - x * sin (2 * π * x) / (4 * π) - cos (2 * π * x) / (8 * π ^ 2)
  have hF : ∀ x ∈ Set.uIcc (0:ℝ) p, HasDerivAt F (x * sin (π * x) ^ 2) x := by
    intro x _
    have h2 : HasDerivAt (fun x : ℝ => 2 * π * x) (2 * π) x := by
      simpa using (hasDerivAt_id x).const_mul (2 * π)
    have hs := h2.sin
    have hc := h2.cos
    have hx2 : HasDerivAt (fun x : ℝ => x ^ 2 / 4) (2 * x / 4) x := by
      simpa using (hasDerivAt_pow 2 x).div_const 4
    have hxs : HasDerivAt (fun x : ℝ => x * sin (2 * π * x) / (4 * π))
        ((1 * sin (2 * π * x) + x * (cos (2 * π * x) * (2 * π))) / (4 * π)) x :=
      ((hasDerivAt_id x).mul hs).div_const (4 * π)
    have hcc : HasDerivAt (fun x : ℝ => cos (2 * π * x) / (8 * π ^ 2))
        ((-sin (2 * π * x) * (2 * π)) / (8 * π ^ 2)) x := hc.div_const (8 * π ^ 2)
    have := (hx2.sub hxs).sub hcc
    convert this using 1
    have hcos : cos (2 * π * x) = 1 - 2 * sin (π * x) ^ 2 := by
      rw [show 2 * π * x = 2 * (π * x) by ring, cos_two_mul, cos_sq']; ring
    rw [hcos]
    field_simp
    ring
  have hint : IntervalIntegrable (fun x : ℝ => x * sin (π * x) ^ 2) MeasureTheory.volume 0 p :=
    (by fun_prop : Continuous (fun x : ℝ => x * sin (π * x) ^ 2)).intervalIntegrable 0 p
  rw [integral_eq_sub_of_hasDerivAt hF hint]
  simp only [F]
  simp
  ring

theorem int_cos_sq (p : ℝ) :
    ∫ x in (0:ℝ)..p, cos (π * x) ^ 2 = p / 2 + sin (2 * π * p) / (4 * π) := by
  have hπ : π ≠ 0 := pi_ne_zero
  let F : ℝ → ℝ := fun x => x / 2 + sin (2 * π * x) / (4 * π)
  have hF : ∀ x ∈ Set.uIcc (0:ℝ) p, HasDerivAt F (cos (π * x) ^ 2) x := by
    intro x _
    have h2 : HasDerivAt (fun x : ℝ => 2 * π * x) (2 * π) x := by
      simpa using (hasDerivAt_id x).const_mul (2 * π)
    have := ((hasDerivAt_id x).div_const 2).add (h2.sin.div_const (4 * π))
    have hcos : cos (2 * π * x) = 2 * cos (π * x) ^ 2 - 1 := by
      rw [show 2 * π * x = 2 * (π * x) by ring, cos_two_mul]
    convert this using 1
    all_goals first
      | (funext y; simp [F])
      | (rw [hcos]; field_simp; ring)
  have hint : IntervalIntegrable (fun x : ℝ => cos (π * x) ^ 2) MeasureTheory.volume 0 p :=
    (by fun_prop : Continuous (fun x : ℝ => cos (π * x) ^ 2)).intervalIntegrable 0 p
  rw [integral_eq_sub_of_hasDerivAt hF hint]
  simp [F]

theorem int_cos (p : ℝ) : ∫ x in (0:ℝ)..p, cos (π * x) = sin (π * p) / π := by
  have hπ : π ≠ 0 := pi_ne_zero
  have := integral_comp_mul_left (fun x => cos x) hπ (a := 0) (b := p)
  simp only [mul_zero] at this
  rw [this, integral_cos]
  simp
  field_simp

end SteklovTriangles
