import Mathlib

/-!
# The fixed-area corollary

`triangle_isoperimetric_heron`: for side lengths `a b c` of a nondegenerate triangle,
`27 * P ≤ (a+b+c)^4` where `P = (a+b+c)(-a+b+c)(a-b+c)(a+b-c)` is Heron's product
(`16 * Area^2 = P`), with equality only for `a = b = c`.  Hence `432 * Area^2 ≤ L^4`,
i.e. `Area / L^2 ≤ √3 / 36`.

`area_corollary`: if `σ * L ≤ S` and `L^2 * k ≥ A` with `k > 0` (here `k = √3/36`),
then `σ * √A ≤ S * √k`; strictness propagates from either hypothesis.
-/

namespace SteklovTriangles

/-- AM–GM for three nonnegative reals, in polynomial form. -/
theorem amgm3 (x y z : ℝ) (hx : 0 ≤ x) (hy : 0 ≤ y) (hz : 0 ≤ z) :
    27 * (x * y * z) ≤ (x + y + z) ^ 3 := by
  nlinarith [sq_nonneg (x - y), sq_nonneg (y - z), sq_nonneg (x - z),
    mul_nonneg hx hy, mul_nonneg hy hz, mul_nonneg hx hz,
    mul_nonneg (mul_nonneg hx hy) hz, sq_nonneg (x + y - 2 * z),
    mul_nonneg hz (sq_nonneg (x - y)), mul_nonneg hx (sq_nonneg (y - z)),
    mul_nonneg hy (sq_nonneg (x - z))]

/-- The isoperimetric inequality for triangles in Heron form. -/
theorem triangle_isoperimetric_heron (a b c : ℝ) (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (h1 : a < b + c) (h2 : b < a + c) (h3 : c < a + b) :
    27 * ((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)) ≤ (a + b + c) ^ 4 := by
  have hx : 0 ≤ -a + b + c := by linarith
  have hy : 0 ≤ a - b + c := by linarith
  have hz : 0 ≤ a + b - c := by linarith
  have key := amgm3 _ _ _ hx hy hz
  have hs : (-a + b + c) + (a - b + c) + (a + b - c) = a + b + c := by ring
  rw [hs] at key
  have hL : 0 ≤ a + b + c := by linarith
  have := mul_le_mul_of_nonneg_left key hL
  nlinarith [this]

/-- Equality in AM–GM (three variables) forces equal arguments. -/
theorem amgm3_eq (x y z : ℝ) (hx : 0 ≤ x) (hy : 0 ≤ y) (hz : 0 ≤ z)
    (h : 27 * (x * y * z) = (x + y + z) ^ 3) : x = y ∧ y = z := by
  -- (x+y+z)^3 - 27xyz ≥ (x+y+z) * ((x-y)^2+(y-z)^2+(z-x)^2) / 2 ≥ 0
  have hsum : 0 ≤ x + y + z := by linarith
  have key : (x + y + z) * ((x - y) ^ 2 + (y - z) ^ 2 + (z - x) ^ 2) ≤ 2 * ((x + y + z) ^ 3 - 27 * (x * y * z)) := by
    nlinarith [sq_nonneg (x - y), sq_nonneg (y - z), sq_nonneg (x - z),
      mul_nonneg hz (sq_nonneg (x - y)), mul_nonneg hx (sq_nonneg (y - z)),
      mul_nonneg hy (sq_nonneg (x - z)), mul_nonneg (mul_nonneg hx hy) hz,
      sq_nonneg (x + y - 2 * z), sq_nonneg (x - 2 * y + z), sq_nonneg (-2 * x + y + z)]
  rw [h, sub_self, mul_zero] at key
  rcases eq_or_lt_of_le hsum with h0 | hpos
  · constructor <;> linarith
  · have hq : (x - y) ^ 2 + (y - z) ^ 2 + (z - x) ^ 2 ≤ 0 := by
      by_contra hcon
      have := mul_pos hpos (lt_of_not_ge hcon)
      linarith
    have e1 : (x - y) ^ 2 = 0 := by nlinarith [sq_nonneg (x - y), sq_nonneg (y - z), sq_nonneg (z - x)]
    have e2 : (y - z) ^ 2 = 0 := by nlinarith [sq_nonneg (x - y), sq_nonneg (y - z), sq_nonneg (z - x)]
    constructor
    · nlinarith [pow_eq_zero_iff (n := 2) (two_ne_zero) |>.mp e1]
    · nlinarith [pow_eq_zero_iff (n := 2) (two_ne_zero) |>.mp e2]

/-- Equality in the Heron-form isoperimetric inequality only for the equilateral triangle. -/
theorem triangle_isoperimetric_heron_eq (a b c : ℝ) (ha : 0 < a) (hb : 0 < b) (hc : 0 < c)
    (h1 : a < b + c) (h2 : b < a + c) (h3 : c < a + b)
    (heq : 27 * ((a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c)) = (a + b + c) ^ 4) :
    a = b ∧ b = c := by
  have hx : 0 ≤ -a + b + c := by linarith
  have hy : 0 ≤ a - b + c := by linarith
  have hz : 0 ≤ a + b - c := by linarith
  have hL : 0 < a + b + c := by linarith
  have hs : (-a + b + c) + (a - b + c) + (a + b - c) = a + b + c := by ring
  have h' : 27 * ((-a + b + c) * (a - b + c) * (a + b - c)) = ((-a + b + c) + (a - b + c) + (a + b - c)) ^ 3 := by
    rw [hs]
    have : (a + b + c) * (27 * ((-a + b + c) * (a - b + c) * (a + b - c))) = (a + b + c) * (a + b + c) ^ 3 := by
      rw [← mul_assoc]; nlinarith [heq]
    exact mul_left_cancel₀ (ne_of_gt hL) this
  obtain ⟨e1, e2⟩ := amgm3_eq _ _ _ hx hy hz h'
  constructor <;> linarith

/-- The fixed-area corollary as real algebra: from `σ L ≤ S` and `A ≤ k L^2`
we get `σ √A ≤ S √k`. -/
theorem area_corollary (σ L A S k : ℝ) (hσ : 0 ≤ σ) (hL : 0 < L) (hk : 0 ≤ k)
    (hper : σ * L ≤ S) (hiso : A ≤ k * L ^ 2) : σ * Real.sqrt A ≤ S * Real.sqrt k := by
  have h1 : Real.sqrt A ≤ Real.sqrt k * L := by
    rw [show Real.sqrt k * L = Real.sqrt (k * L ^ 2) by
      rw [Real.sqrt_mul hk, Real.sqrt_sq hL.le]]
    exact Real.sqrt_le_sqrt hiso
  calc σ * Real.sqrt A ≤ σ * (Real.sqrt k * L) := mul_le_mul_of_nonneg_left h1 hσ
    _ = (σ * L) * Real.sqrt k := by ring
    _ ≤ S * Real.sqrt k := mul_le_mul_of_nonneg_right hper (Real.sqrt_nonneg k)

/-- Strict version: strictness of the perimeter inequality gives a strict area inequality. -/
theorem area_corollary_strict (σ L A S k : ℝ) (hσ : 0 ≤ σ) (hL : 0 < L) (hk : 0 < k)
    (hper : σ * L < S) (hiso : A ≤ k * L ^ 2) : σ * Real.sqrt A < S * Real.sqrt k := by
  have h1 : Real.sqrt A ≤ Real.sqrt k * L := by
    rw [show Real.sqrt k * L = Real.sqrt (k * L ^ 2) by
      rw [Real.sqrt_mul hk.le, Real.sqrt_sq hL.le]]
    exact Real.sqrt_le_sqrt hiso
  have hsk : 0 < Real.sqrt k := Real.sqrt_pos.mpr hk
  calc σ * Real.sqrt A ≤ σ * (Real.sqrt k * L) := mul_le_mul_of_nonneg_left h1 hσ
    _ = (σ * L) * Real.sqrt k := by ring
    _ < S * Real.sqrt k := mul_lt_mul_of_pos_right hper hsk

end SteklovTriangles
