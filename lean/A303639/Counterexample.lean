/-
Copyright 2026 Sunsu Jeong (github.com/DCLXAI).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-/
import Mathlib.NumberTheory.SumTwoSquares
import Mathlib.Tactic

/-!
# A counterexample to Zhi-Wei Sun's conjecture on OEIS A303639

OEIS [A303639](https://oeis.org/A303639) counts the representations of `n` as
`a ^ 2 + b ^ 2 + C(2 * c + 1, c) + C(2 * d + 1, d)` with `a, b, c, d` nonnegative integers.
Zhi-Wei Sun conjectured that this count is positive for every `n > 1`, and verified it
for `n` up to `6 * 10 ^ 8`.

This file proves the conjecture is **false**: no such representation exists for
`n = 800322180`. Writing `B k = C(2 * k + 1, k)`, the proof runs as follows.

* `B` is monotone and `B 16 = 1166803110 > 800322180`, so only `c, d ≤ 15` can occur
  (`A303639.lt_sixteen`).
* That leaves the 136 remainders `r = 800322180 - B c - B d`, listed in `A303639.certs`
  together with a prime `p ≡ 3 [MOD 4]` dividing `r` to an odd exponent `e`.
* By Fermat's two-square theorem (`Nat.eq_sq_add_sq_iff`) such an `r` is not a sum of two
  squares (`A303639.not_sq_add_sq`), so no representation exists (`A303639.no_representation`).

The counterexample was found by an exhaustive search and independently reconfirmed by three
implementations; it is recorded as an approved comment on the OEIS entry.

*References:*
- [A303639](https://oeis.org/A303639)
- Search, certificates and independent verifiers:
  <https://github.com/DCLXAI/a303639-counterexample> (DOI: 10.5281/zenodo.21863025)
-/

namespace A303639

/-- `B k` is the binomial term `C(2 * k + 1, k)` appearing in A303639. -/
def B (k : ℕ) : ℕ := (2 * k + 1).choose k

theorem B_le_succ (k : ℕ) : B k ≤ B (k + 1) := by
  have e1 : B (k + 1) = (2 * k + 2).choose k + (2 * k + 2).choose (k + 1) := by
    show (2 * (k + 1) + 1).choose (k + 1) = _
    have h : 2 * (k + 1) + 1 = (2 * k + 2) + 1 := by ring
    rw [h, Nat.choose_succ_succ]
  have e2 : (2 * k + 2).choose (k + 1) = (2 * k + 1).choose k + (2 * k + 1).choose (k + 1) := by
    have h : 2 * k + 2 = (2 * k + 1) + 1 := by ring
    rw [h, Nat.choose_succ_succ]
  simp only [B] at e1 e2 ⊢
  omega

theorem B_mono : Monotone B := monotone_nat_of_le_succ B_le_succ

theorem B_0 : B 0 = 1 := by
  show (1).choose 0 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_1 : B 1 = 3 := by
  show (3).choose 1 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_2 : B 2 = 10 := by
  show (5).choose 2 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_3 : B 3 = 35 := by
  show (7).choose 3 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_4 : B 4 = 126 := by
  show (9).choose 4 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_5 : B 5 = 462 := by
  show (11).choose 5 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_6 : B 6 = 1716 := by
  show (13).choose 6 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_7 : B 7 = 6435 := by
  show (15).choose 7 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_8 : B 8 = 24310 := by
  show (17).choose 8 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_9 : B 9 = 92378 := by
  show (19).choose 9 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_10 : B 10 = 352716 := by
  show (21).choose 10 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_11 : B 11 = 1352078 := by
  show (23).choose 11 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_12 : B 12 = 5200300 := by
  show (25).choose 12 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_13 : B 13 = 20058300 := by
  show (27).choose 13 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_14 : B 14 = 77558760 := by
  show (29).choose 14 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_15 : B 15 = 300540195 := by
  show (31).choose 15 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_16 : B 16 = 1166803110 := by
  show (33).choose 16 = _
  rw [Nat.choose_eq_descFactorial_div_factorial]
  norm_num [Nat.descFactorial, Nat.factorial]

theorem B_ge (k : ℕ) (hk : 16 ≤ k) : 1166803110 ≤ B k := by
  have := B_mono hk
  rwa [B_16] at this

/-- Since `B 16 > 800322180`, any representation of `800322180` has `c, d ≤ 15`. -/
theorem lt_sixteen {a b c d : ℕ} (h : 800322180 = a ^ 2 + b ^ 2 + B c + B d) :
    c ≤ 15 ∧ d ≤ 15 := by
  constructor <;> by_contra hk <;>
    [have := B_ge c (by omega); have := B_ge d (by omega)] <;> omega

/--
Fermat's two-square theorem in certificate form: if some prime `p ≡ 3 [MOD 4]` divides `r`
to an odd exponent `e`, then `r` is not a sum of two squares.
-/
theorem not_sq_add_sq (r p e : ℕ) (hr : r ≠ 0) (hp : p.Prime) (hp4 : p % 4 = 3)
    (hdvd : p ^ e ∣ r) (hnd : ¬ p ^ (e + 1) ∣ r) (he : Odd e) :
    ¬ ∃ x y : ℕ, r = x ^ 2 + y ^ 2 := by
  rw [Nat.eq_sq_add_sq_iff]
  push_neg
  have hfac : r.factorization p = e := by
    have h1 : e ≤ r.factorization p := (hp.pow_dvd_iff_le_factorization hr).1 hdvd
    have h2 : ¬ (e + 1 ≤ r.factorization p) := fun h =>
      hnd ((hp.pow_dvd_iff_le_factorization hr).2 h)
    omega
  have hpdvd : p ∣ r := dvd_trans (dvd_pow_self p (by have := he.pos; omega)) hdvd
  refine ⟨p, Nat.mem_primeFactors.2 ⟨hp, hpdvd, hr⟩, hp4, ?_⟩
  rw [← Nat.factorization_def _ hp, hfac]
  simpa [Nat.even_iff, Nat.odd_iff] using he

/--
The 136 certificates `(r, p, e)`: `r` runs over `800322180 - B c - B d` for `0 ≤ c ≤ d ≤ 15`,
and `p ^ e` exactly divides `r` with `p ≡ 3 [MOD 4]` prime and `e` odd.
-/
def certs : List (ℕ × ℕ × ℕ) :=
  [(800322178, 647, 1),
   (800322176, 163, 1),
   (800322169, 67, 1),
   (800322144, 3, 1),
   (800322053, 23, 1),
   (800321717, 467, 1),
   (800320463, 223, 1),
   (800315744, 991, 1),
   (800297869, 7, 1),
   (800229801, 3, 1),
   (799969463, 23, 1),
   (798970101, 3, 1),
   (795121879, 887, 1),
   (780263879, 26905651, 1),
   (722763419, 7, 1),
   (499781984, 1999, 1),
   (800322174, 3, 3),
   (800322167, 643, 1),
   (800322142, 19, 1),
   (800322051, 3, 1),
   (800321715, 3, 7),
   (800320461, 3, 1),
   (800315742, 3, 1),
   (800297867, 78607, 1),
   (800229799, 27594131, 1),
   (799969461, 3, 1),
   (798970099, 31, 1),
   (795121877, 11, 1),
   (780263877, 3, 1),
   (722763417, 3, 1),
   (499781982, 3, 1),
   (800322160, 11, 1),
   (800322135, 3, 1),
   (800322044, 200080511, 1),
   (800321708, 23, 1),
   (800320454, 284003, 1),
   (800315735, 4326031, 1),
   (800297860, 19, 1),
   (800229792, 3, 1),
   (799969454, 399984727, 1),
   (798970092, 3, 1),
   (795121870, 79512187, 1),
   (780263870, 127, 1),
   (722763410, 2819, 1),
   (499781975, 7, 1),
   (800322110, 7, 1),
   (800322019, 7, 1),
   (800321683, 7, 1),
   (800320429, 43, 1),
   (800315710, 5843, 1),
   (800297835, 3, 1),
   (800229767, 7923067, 1),
   (799969429, 7, 1),
   (798970067, 7, 1),
   (795121845, 3, 1),
   (780263845, 23, 1),
   (722763385, 23, 1),
   (499781950, 23, 1),
   (800321928, 3, 1),
   (800321592, 3, 1),
   (800320338, 719, 1),
   (800315619, 3, 1),
   (800297744, 563, 1),
   (800229676, 200057419, 1),
   (799969338, 3, 3),
   (798969976, 7, 1),
   (795121754, 7, 1),
   (780263754, 3, 1),
   (722763294, 3, 1),
   (499781859, 3, 1),
   (800321256, 3, 3),
   (800320002, 3, 1),
   (800315283, 3, 1),
   (800297408, 12504647, 1),
   (800229340, 23, 1),
   (799969002, 3, 1),
   (798969640, 7, 1),
   (795121418, 7, 1),
   (780263418, 3, 1),
   (722762958, 3, 1),
   (499781523, 3, 1),
   (800318748, 3, 1),
   (800314029, 19, 1),
   (800296154, 7, 1),
   (800228086, 7, 1),
   (799967748, 3, 1),
   (798968386, 139, 1),
   (795120164, 227, 1),
   (780262164, 11, 1),
   (722761704, 3, 5),
   (499780269, 67, 1),
   (800309310, 3, 1),
   (800291435, 199, 1),
   (800223367, 31, 1),
   (799963029, 31, 1),
   (798963667, 3331, 1),
   (795115445, 139, 1),
   (780257445, 3, 1),
   (722756985, 3, 1),
   (499775550, 3, 1),
   (800273560, 689891, 1),
   (800205492, 3, 1),
   (799945154, 223, 1),
   (798945792, 3, 1),
   (795097570, 59, 1),
   (780239570, 11, 1),
   (722739110, 7591, 1),
   (499757675, 19990307, 1),
   (800137424, 19, 1),
   (799877086, 887, 1),
   (798877724, 19, 1),
   (795029502, 3, 1),
   (780171502, 11, 1),
   (722671042, 19, 1),
   (499689607, 19, 1),
   (799616748, 3, 1),
   (798617386, 19, 1),
   (794769164, 7, 1),
   (779911164, 19, 1),
   (722410704, 3, 3),
   (499429269, 19, 1),
   (797618024, 7, 1),
   (793769802, 3, 1),
   (778911802, 19, 1),
   (721411342, 19, 1),
   (498429907, 19, 1),
   (789921580, 7, 1),
   (775063580, 19, 1),
   (717563120, 19, 1),
   (494581685, 19, 1),
   (760205580, 3, 1),
   (702705120, 3, 1),
   (479723685, 3, 1),
   (645204660, 3, 1),
   (422223225, 3, 1),
   (199241790, 3, 1)]

theorem certs_ok : ∀ t ∈ certs, ¬ ∃ x y : ℕ, t.1 = x ^ 2 + y ^ 2 := by
  rintro ⟨r, p, e⟩ ht
  refine not_sq_add_sq r p e ?_ ?_ ?_ ?_ ?_ ?_ <;>
    · fin_cases ht <;> norm_num

theorem not_sq_add_sq_of_mem (r : ℕ) (hr : r ∈ certs.map Prod.fst) :
    ¬ ∃ x y : ℕ, r = x ^ 2 + y ^ 2 := by
  obtain ⟨t, ht, rfl⟩ := List.mem_map.1 hr
  exact certs_ok t ht

theorem mem_certs (c d : ℕ) (hc : c ≤ 15) (hd : d ≤ 15) :
    800322180 - B c - B d ∈ certs.map Prod.fst := by
  interval_cases c <;> interval_cases d <;>
    simp only [B_0, B_1, B_2, B_3, B_4, B_5, B_6, B_7,
      B_8, B_9, B_10, B_11, B_12, B_13, B_14, B_15] <;>
    decide

/--
**A303639 is zero at `n = 800322180`**: there is no way to write `800322180` as
`a ^ 2 + b ^ 2 + C(2 * c + 1, c) + C(2 * d + 1, d)` with `a, b, c, d` nonnegative integers.
-/
theorem no_representation : ¬ ∃ a b c d : ℕ,
    800322180 = a ^ 2 + b ^ 2 + (2 * c + 1).choose c + (2 * d + 1).choose d := by
  rintro ⟨a, b, c, d, h⟩
  rw [show (2 * c + 1).choose c = B c from rfl, show (2 * d + 1).choose d = B d from rfl] at h
  obtain ⟨hc, hd⟩ := lt_sixteen h
  exact not_sq_add_sq_of_mem _ (mem_certs c d hc hd) ⟨a, b, by omega⟩

/--
**Sun's conjecture on A303639 is false.** Sun conjectured that every `n > 1` is representable as
`a ^ 2 + b ^ 2 + C(2 * c + 1, c) + C(2 * d + 1, d)`; `n = 800322180` is a counterexample.
-/
theorem sun_conjecture_false :
    ¬ ∀ n : ℕ, 1 < n → ∃ a b c d : ℕ,
      n = a ^ 2 + b ^ 2 + (2 * c + 1).choose c + (2 * d + 1).choose d := by
  intro h
  exact no_representation (h 800322180 (by norm_num))

end A303639
