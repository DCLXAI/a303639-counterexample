# Lean 4 formalisation

A machine-checked proof, in Lean 4 with mathlib, that Zhi-Wei Sun's conjecture on
[OEIS A303639](https://oeis.org/A303639) is false.

## Statement

`A303639/Counterexample.lean` proves, with no `sorry` and no `native_decide`:

```lean
theorem no_representation :
    ¬ ∃ a b c d : ℕ, 800322180 = a ^ 2 + b ^ 2 + (2 * c + 1).choose c + (2 * d + 1).choose d

theorem sun_conjecture_false :
    ¬ ∀ n : ℕ, 1 < n → ∃ a b c d : ℕ,
      n = a ^ 2 + b ^ 2 + (2 * c + 1).choose c + (2 * d + 1).choose d
```

The statement drops Sun's normalisations `a ≤ b` and `c ≤ d`. Those do not change whether a
representation exists, and refuting the unordered form is the stronger result: a representation
with `a ≤ b`, `c ≤ d` is in particular a representation.

## Proof outline

Write `B k = C(2k+1, k)`.

1. `B_mono` — `B` is monotone, by Pascal's rule (`Nat.choose_succ_succ`).
2. `lt_sixteen` — `B 16 = 1166803110 > 800322180`, so any representation has `c, d ≤ 15`.
3. `certs` — the resulting 136 remainders `r = 800322180 - B c - B d` (`0 ≤ c ≤ d ≤ 15`), each
   paired with a prime `p ≡ 3 (mod 4)` and an odd exponent `e` with `p^e ‖ r`. These are the
   certificates in [`../certificates/certificates.json`](../certificates/certificates.json),
   with `p` chosen as the smallest valid certifying prime for each remainder.
4. `not_sq_add_sq` — Fermat's two-square theorem (mathlib's `Nat.eq_sq_add_sq_iff`) turns each
   certificate into a proof that `r` is not a sum of two squares.
5. `no_representation` — combining the above, no representation of `800322180` exists.

## Building

```sh
cd lean
lake exe cache get   # downloads prebuilt mathlib oleans
lake build
```

Pinned to Lean `v4.27.0` and mathlib `v4.27.0` (rev `a3a10db0e9d66acbebf76c5e6a135066525ac900`),
matching the toolchain of
[google-deepmind/formal-conjectures](https://github.com/google-deepmind/formal-conjectures),
where the corresponding statement lives as `FormalConjectures/OEIS/303639.lean`.

To confirm that nothing is assumed beyond Lean's standard axioms:

```lean
#print axioms A303639.sun_conjecture_false
-- 'A303639.sun_conjecture_false' depends on axioms: [propext, Classical.choice, Quot.sound]
```
