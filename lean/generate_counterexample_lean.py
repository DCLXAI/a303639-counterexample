#!/usr/bin/env python3
"""Regenerate `A303639/Counterexample.lean`.

Recomputes the 136 non-sum-of-two-squares certificates for n = 800322180 from scratch,
cross-checks them against ../certificates/certificates.json, and emits the Lean 4 proof.

Usage:  python3 generate_counterexample_lean.py > A303639/Counterexample.lean
"""
import json, os, sys
from math import comb

N = 800322180
B = lambda k: comb(2 * k + 1, k)


def factor(m):
    f = {}
    d = 2
    while d * d <= m:
        while m % d == 0:
            f[d] = f.get(d, 0) + 1
            m //= d
        d += 1 if d == 2 else 2
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return f


def is_prime(x):
    if x < 2:
        return False
    if x % 2 == 0:
        return x == 2
    d = 3
    while d * d <= x:
        if x % d == 0:
            return False
        d += 2
    return True


# smallest certifying prime p = 3 (mod 4) with odd multiplicity, for every pair 0 <= c <= d <= 15
items = []
for c in range(16):
    for d in range(c, 16):
        r = N - B(c) - B(d)
        assert r > 0
        f = factor(r)
        cand = sorted(p for p, e in f.items() if p % 4 == 3 and e % 2 == 1)
        assert cand, (c, d, r)
        p, e = cand[0], f[cand[0]]
        assert is_prime(p) and r % p ** e == 0 and r % p ** (e + 1) != 0
        items.append((c, d, r, p, e))
assert len(items) == 136

# cross-check against the published certificates
pub_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "certificates", "certificates.json")
pub = json.load(open(pub_path))
assert pub["n"] == N and pub["count"] == 136
by_pair = {(c, d): r for (c, d, r, _, _) in items}
for x in pub["certificates"]:
    assert by_pair[(x["c"], x["d"])] == x["r"], x
    assert is_prime(x["p"]) and x["p"] % 4 == 3 and x["e"] % 2 == 1
    assert x["r"] % x["p"] ** x["e"] == 0 and x["r"] % x["p"] ** (x["e"] + 1) != 0
print("cross-check against certificates.json: OK (136 pairs)", file=sys.stderr)

bvals = "\n".join(
    f"theorem B_{k} : B {k} = {B(k)} := by\n"
    f"  show ({2*k+1}).choose {k} = _\n"
    f"  rw [Nat.choose_eq_descFactorial_div_factorial]\n"
    f"  norm_num [Nat.descFactorial, Nat.factorial]\n"
    for k in range(17))

certlist = ",\n   ".join(f"({r}, {p}, {e})" for (_, _, r, p, e) in items)
bsimp = ", ".join(f"B_{k}" for k in range(8)) + ",\n      " + ", ".join(f"B_{k}" for k in range(8, 16))

print(f"""/-
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

{bvals}
theorem B_ge (k : ℕ) (hk : 16 ≤ k) : 1166803110 ≤ B k := by
  have := B_mono hk
  rwa [B_16] at this

/-- Since `B 16 > 800322180`, any representation of `800322180` has `c, d ≤ 15`. -/
theorem lt_sixteen {{a b c d : ℕ}} (h : 800322180 = a ^ 2 + b ^ 2 + B c + B d) :
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
  [{certlist}]

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
    simp only [{bsimp}] <;>
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

end A303639""")
