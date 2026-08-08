# A counterexample to a representation conjecture of Sun (OEIS A303639)
*(research note skeleton — 4-6쪽 분량 목표)*

## 1. The conjecture
OEIS A303639 (Zhi-Wei Sun, Apr 27 2018) counts representations
a(n) = #{ (a,b,c,d) : n = a^2 + b^2 + C(2c+1,c) + C(2d+1,d), a <= b, c <= d, all >= 0 }.
The entry states: "Conjecture: a(n) > 0 for all n > 1", noting it "is similar to
the author's conjecture in A303540" and that "it has been verified that a(n) > 0
for all n = 2..6*10^8". Related work: Z.-W. Sun, "Refining Lagrange's four-square
theorem", J. Number Theory 175 (2017) 167-190; "Restricted sums of four squares",
arXiv:1701.05868.

## 2. Result
n = 800322180 is the unique exception in (1, 2*10^9]. Certificate: 136 pairs, each
remainder fails Fermat's two-squares criterion; machine-checkable JSON + stdlib checker.

## 3. Method and performance
Bitset shift-OR coverage over [1, N], N = 2*10^9, single core:
base = { a^2 + b^2 <= N } filled by the direct double loop (~pi/8 * N ~ 7.9*10^8
iterations, 15.2 s); coverage = OR of 136 copies of base shifted by B(c)+B(d)
(8.0 s); linear scan for uncovered n (1.2 s). Total 24.4 s, 500 MB (two bitsets),
gcc -O3. Output: the single uncovered value n = 800322180.
The result has three independent implementations: the original discovery search,
a from-scratch rewrite (searcher/search.c, timings above), and the
factorization-based verifier (verifier/verify.py); all agree.

## 4. Why the counterexample exists here — and why it is probably unique

### 4.1 A rigorous local obstruction (Lemma)

**Lemma 4.1.** Let M = 8*9*49*361 = 1,273,608 and suppose n == 800322180 (mod M).
Let B(k) = binomial(2k+1,k) and let (c,d), c <= d, be any of the 106 pairs listed
in `certificates/forced_table.json` with B(c)+B(d) <= n. Then
r = n - B(c) - B(d) is not a sum of two squares. Consequently, for such n,
a representation n = a^2 + b^2 + B(c) + B(d) can only come from the remaining
30 pairs.

*Proof.* r mod M is determined by n mod M, and the table lists r mod M for each
pair. Each entry falls into one of four cases, each of which forbids r from being
a sum of two squares by Fermat's theorem (m is a sum of two squares iff every
prime p == 3 (mod 4) divides m to even order):

- **(mod8, 51 pairs)** r == 3, 6 or 7 (mod 8). Then the odd part m of r satisfies
  m == 3 (mod 4), so the number of primes == 3 (mod 4) dividing m counted with
  multiplicity is odd; some such prime occurs to odd order.
- **(v3, 29 pairs)** r == 3 or 6 (mod 9): then 3 | r and 9 !| r, so v_3(r) = 1.
- **(v7, 13 pairs)** r == 7u (mod 49) with u a unit: v_7(r) = 1.
- **(v19, 13 pairs)** r == 19u (mod 361) with u a unit: v_19(r) = 1.

The classification of all 106 pairs is a finite residue computation, generated
and machine-verified by `analysis/forced_table.py`. QED

The counterexample n0 = 800322180 = 2^2 * 3 * 5 * 7 * 19 * 100291 satisfies the
hypothesis, so only its 30 free pairs could have produced a representation; the
certificates show none did.

### 4.2 How extreme is this class? (computation)

Over 3*10^5 uniformly random n in [4*10^8, 2*10^9], the number F(n) of
non-forced ("free") pairs averages 55.1 with observed minimum exactly 30:
n0's class sits at the 0.005th percentile of obstruction — essentially the most
obstructed residue class this problem admits at this size.

### 4.3 Heuristic probability model (not rigorous)

Marginal Landau-Ramanujan density of sums of two squares near 8*10^8 is
~0.7642/sqrt(ln x) ~ 0.17 (measured 0.174); naive independence gives
P(all 136 fail) ~ 4.9e-12 and E[#exceptions in [4e8, 2e9]] ~ 0.025.
Conditioning on n0's class: the 106 forced pairs fail with probability 1
(Lemma 4.1), and the free-pair conditional density rises to d_cond ~ 0.428,
giving P(exception | class) ~ 3.6e-8 — an amplification of ~7*10^3.
Validation: for 300 random n in the same class the model predicts a mean of
122.3 failing pairs; observed 122.4 (max 131/136; no full exception, as
expected). Refined E[#exceptions in [4e8, 2e9]] ~ 0.053; observing exactly one
is a ~5% event.

### 4.4 Uniqueness (heuristic conjecture)

The per-n exception probability decays superpolynomially — K(n) ~ (log_4 n)^2/2
qualifying pairs against a density decaying only as 1/sqrt(ln n) — so the tail
expectation converges: E[#exceptions, n > 2*10^9] ~ 0.007 (naive) to ~0.014
(refined), essentially exhausted by 10^12.

**Conjecture.** n = 800322180 is the only integer n > 1 with a(n) = 0.

**Corollary (targeted search).** Any search beyond 2*10^9 need only scan classes
with small F(n), a >100x reduction of the search space.

*(Caveats: 4.3-4.4 rest on independence heuristics beyond the forced classes;
class-probability sampling error ~2x; all numbers reproducible via
analysis/model.py and analysis/forced_table.py, seeded.)*

## 5. Verification and reproducibility
Repo layout, certificate format, checker guarantees (completeness + validity).
Two further consistency checks against the OEIS entry itself: (i) our
implementation reproduces the official DATA (all 80 published terms) and every
a(n) = 1 EXAMPLE value (n = 2530, 3258, 5300, 13453, 20964); (ii) the
certificate criterion used here (a prime p == 3 (mod 4) to odd order) is the
same sum-of-two-squares test as the QQ[] function in the entry's own Mathematica
program, so the author's and our verification logic coincide.

## 6. Remarks
Relation to Crocker/Sun-type mixed representation problems. AI-assisted workflow note.
