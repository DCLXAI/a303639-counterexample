# A counterexample to a representation conjecture of Sun (OEIS A303639)
*(research note — 4-6쪽 분량 목표)*

## 1. The conjecture
OEIS A303639 (Zhi-Wei Sun, Apr 27 2018) counts representations
a(n) = #{ (a,b,c,d) : n = a^2 + b^2 + C(2c+1,c) + C(2d+1,d), a <= b, c <= d, all >= 0 }.
The entry states: "Conjecture: a(n) > 0 for all n > 1", noting it "is similar to
the author's conjecture in A303540" and that "it has been verified that a(n) > 0
for all n = 2..6*10^8". Related work: Z.-W. Sun, "Refining Lagrange's four-square
theorem", J. Number Theory 175 (2017) 167-190; "Restricted sums of four squares",
arXiv:1701.05868.

## 2. Result

Write B(k) = C(2k+1,k).

**Theorem 1.** a(800322180) = 0. Equivalently, n0 = 800322180 is not of the form
a^2 + b^2 + B(c) + B(d). This refutes the conjecture.

*Proof.* B(16) = 1166803110 > n0, so any representation has c <= d <= 15; there
are exactly 136 pairs (c,d) with 0 <= c <= d <= 15, and all 136 satisfy
B(c) + B(d) <= n0. For each, `certificates/certificates.json` supplies a prime
p == 3 (mod 4) and an odd e with p^e || r, r = n0 - B(c) - B(d). By Fermat's
two-squares theorem no such r is a sum of two squares, so no pair admits a
representation. QED

The proof is finite, explicit and machine-checkable: `check_certificates.py`
re-derives the pair list from n0 and validates all 136 certificates using only
the Python standard library (Section 5).

**Proposition 2 (computational).** n0 is the only n with 1 < n <= 2*10^9 and
a(n) = 0.

Proposition 2 rests on the exhaustive scan of Section 3, not on certificates,
and is stated as a computational result.

## 3. Method and performance
Bitset shift-OR coverage over [1, N], N = 2*10^9, single core:
base = { a^2 + b^2 <= N } filled by the direct double loop (~pi/8 * N ~ 7.9*10^8
iterations); coverage = OR of one shifted copy of base per qualifying pair
(shift B(c)+B(d)); linear scan for uncovered n.

The number of qualifying pairs at N = 2*10^9 is **152**, not 136: the pair count
is governed by B(k) <= N, and B(16) = 1166803110 is itself below 2*10^9, so
index 16 pairs with every index 0..15. (Only (16,16) is excluded, since
2*B(16) > N.) The 136 figure is the count relevant to n0 alone, where the
binding constraint is B(16) > n0. `searcher/search.c` now generates B(k) while
B(k) <= N rather than hard-coding a cutoff; see the erratum in Section 5.6.

Measured on an Apple M4, clang -O3, single core: base fill 4.2 s, 152 shifts
1.3 s, scan 0.6 s — 6.1 s total, 500 MB (two bitsets). Output: the single
uncovered value n = 800322180.

The result has three independent implementations: the original discovery search,
a from-scratch rewrite (`searcher/search.c`, timings above), and the
factorization-based verifier (`verifier/verify.py`); all agree.

## 4. Why the counterexample exists here — and why it is probably unique

### 4.1 A rigorous local obstruction (Lemma)

**Lemma 4.1.** Let M = 8*9*49*361 = 1,273,608 and suppose n == 800322180 (mod M).
Let (c,d), c <= d, be any of the 106 pairs listed in
`certificates/forced_table.json`. Then r = n - B(c) - B(d) is not a sum of two
squares. If moreover n < B(16) + B(0) = 1,166,803,111, those 106 pairs together
with the 30 listed as free exhaust the pairs with B(c)+B(d) <= n, so a
representation n = a^2 + b^2 + B(c) + B(d) can only come from the 30 free pairs.

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

The first assertion depends only on r mod M and so holds for every n in the
class. The second is what needs the bound n < B(16)+B(0): past that point the
pair (0,16) qualifies and the 136-pair list is no longer complete, so "only the
30 free pairs remain" would be false. (This bound was omitted in v1.0.0; see
Section 5.6.)

The counterexample n0 = 800322180 = 2^2 * 3 * 5 * 7 * 19 * 100291 satisfies both
hypotheses — n0 < 1,166,803,111 — so only its 30 free pairs could have produced a
representation; the certificates show none did.

### 4.2 How extreme is this class? (computation)

Over 3*10^5 uniformly random n in [4*10^8, 2*10^9], the number F(n) of
non-forced ("free") pairs averages **58.32** with observed minimum exactly 30:
n0's class sits at the **0.0023rd percentile** of obstruction — essentially the
most obstructed residue class this problem admits at this size. Over the same
sample the total qualifying-pair count K(n) averages 143.9 and ranges over
135..152, the jump reflecting the index-16 pairs entering above 1,166,803,111.

### 4.3 Heuristic probability model (not rigorous)

Marginal Landau-Ramanujan density of sums of two squares near 8*10^8 is
~0.7642/sqrt(ln x) ~ 0.17 (measured 0.174); naive independence over n0's 136
pairs gives P(all fail) ~ 5.1e-12, and averaging (1 - 0.174)^K(n) over the
sample gives E[#exceptions in [4e8, 2e9]] ~ 0.004.
Conditioning on n0's class: the 106 forced pairs fail with probability 1
(Lemma 4.1), and the free-pair conditional density rises to d_cond ~ 0.424,
giving P(exception | class) ~ 6.5e-8 — an amplification of ~1.3*10^4.
Validation: for 300 random n in the same class the model predicts a mean of
128.6 failing pairs; observed 129.0 out of a mean 142.0 qualifying pairs
(max 143; no full exception, as expected). Refined
E[#exceptions in [4e8, 2e9]] ~ 0.041; observing exactly one is a ~4% event.

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
analysis/model.py and analysis/forced_table.py, seeded. Nothing in 4.2-4.4
enters the proof of Theorem 1.)*

## 5. Verification and reproducibility

### 5.1 Artifact layout
- `searcher/search.c` — exhaustive bitset scan of [1, N] (Proposition 2).
- `verifier/verify.py` — independent check of Theorem 1 by direct factorization
  (sympy); seconds.
- `certificates/certificates.json` — the 136 certificates.
- `certificates/check_certificates.py` — trustless checker, standard library only.
- `analysis/forced_table.py` — generates and verifies the Lemma 4.1 table.
- `analysis/obstruction.py`, `analysis/model.py` — Section 4 statistics.
- `tests/` — pytest suite; `.github/workflows/test.yml` — CI.

### 5.2 Certificate format
`certificates.json` holds `n`, `count`, and a list of records
`{c, d, Bc, Bd, r, p, e}`: the pair indices, the two binomial values, the
remainder r = n - Bc - Bd, and the certifying prime power p^e. Example:
`{"c": 0, "d": 0, "Bc": 1, "Bd": 1, "r": 800322178, "p": 647, "e": 1}`.
Certifying primes are small — the largest is 399,984,727 — so every field is
checkable by hand in principle.

### 5.3 What the checker guarantees
`check_certificates.py` establishes both halves of the proof:

*Completeness.* It does not trust the supplied pair list. It regenerates
B(0), B(1), ... from `math.comb` until B(k) > n, takes KMAX as the last index
with B(KMAX) <= n (asserting B(KMAX) <= n < B(KMAX+1), which yields KMAX = 15
and B(16) = 1166803110 for this n), enumerates every (i,j) with
i <= j <= KMAX and B(i)+B(j) <= n, and asserts that this set equals the set of
certificate keys exactly — no pair missing, no extra pair.

*Validity.* For each pair it recomputes Bc, Bd and r from n, then checks
p == 3 (mod 4), e odd, p prime, p^e | r and p^(e+1) !| r.

Primality is by unconditional trial division to isqrt(p), not Miller-Rabin.
This is deliberate: a Miller-Rabin checker would make correctness depend on a
published bound for the smallest strong pseudoprime to the chosen base set,
which a referee must then look up and audit. (v1.0.0 cited such a bound
incorrectly; Section 5.6.) Trial division needs no external input, and since
p < 2^31 the entire run takes about 0.02 s.

No factorization is performed anywhere in the checker, and it imports only
`json`, `sys` and `math.comb`.

### 5.4 Independent implementations
Theorem 1 is confirmed by three routes that share no code: the certificate
checker above (residue arithmetic only), `verifier/verify.py` (full sympy
factorization of each r), and the searcher, which rediscovers n0 without being
told about it. Proposition 2 has been reproduced by the original discovery
search and by the from-scratch `search.c`.

### 5.5 Cross-checks against the OEIS entry itself
Two checks guard against having implemented a different sequence than A303639:

1. A direct-definition implementation reproduces the entry's published DATA —
   all 80 terms — and every a(n) = 1 value listed in its EXAMPLE section
   (n = 2530, 3258, 5300, 13453, 20964). Both are asserted in `tests/`.
2. The certificate criterion used here (a prime p == 3 (mod 4) to odd order) is
   the same sum-of-two-squares test as the `QQ[]` function in the entry's own
   Mathematica program, so the author's verification logic and ours coincide.

### 5.6 Erratum for v1.0.0
Release v1.0.0 (Zenodo DOI 10.5281/zenodo.21863025) contains four defects, all
now corrected. None affects Theorem 1, whose certificates were valid as
published; the searcher and the Section 4 statistics were affected.

1. **Searcher pair cutoff.** `search.c` hard-coded indices 0..15 on the strength
   of the false comment "B(15) <= 2e9 < B(16)". In fact B(16) = 1166803110 <
   2*10^9, so the scan used 136 shifts where 152 are required. Proposition 2
   nevertheless survives unchanged: the restricted scan covers a subset of the
   true coverage set, so the n it certifies as covered are genuinely covered,
   and the one value it left uncovered, n0, is below B(16) and therefore could
   never have used an index-16 pair. The corrected 152-shift searcher has been
   run over the full range and returns the same unique value. The searcher now
   derives its cutoff from N, and a regression test pins B(16), the 152-pair
   count, and the searcher's behaviour at N = B(16)+B(0).
2. **Statistics pair list.** `analysis/model.py` fixed its pair list at n0's 136
   pairs while sampling n up to 2*10^9, so it measured a restricted problem
   above 1,166,803,111. All Section 4.2-4.4 figures have been recomputed with
   the pair list derived per n; the revised values appear above. The
   qualitative conclusion — n0 sits in an extraordinarily obstructed class, and
   the expected number of further exceptions is small — is unchanged.
3. **Lemma 4.1 scope.** The lemma's conclusion was stated for all n == n0
   (mod M) without the bound n < B(16)+B(0) that it requires. Added.
4. **Primality comment.** The checker's Miller-Rabin bases 2..37 were annotated
   "deterministic for m < 3.3*10^24". That bound belongs to the 13 bases 2..41
   (smallest strong pseudoprime 3317044064679887385961981); for the 12 bases
   actually used it is 318665857834031151167461 ~ 3.19*10^23
   (Sorenson-Webster, "Strong pseudoprimes to twelve prime bases"). It was
   never binding, since p < 2^31, but the test has been replaced by trial
   division so that no such bound is needed.

Anyone citing this work should cite the corrected release rather than v1.0.0.

### 5.7 Reproducing
```
python3 certificates/check_certificates.py certificates/certificates.json
python3 verifier/verify.py          # requires sympy
python3 analysis/forced_table.py
python3 analysis/model.py           # requires numpy + sympy; ~2-3 min
pytest tests/ -q
gcc -O3 -o search searcher/search.c && ./search 2000000000
```
CI runs all of these on every push, including the full 2*10^9 scan.

## 6. Remarks

### 6.1 Context
A303639 belongs to a family of conjectures of Sun asserting that every
sufficiently large integer is a sum of two squares plus two terms drawn from a
sparse sequence — here the central-binomial-type values B(k) = C(2k+1,k). Such
statements are plausible on density grounds: the sums of two squares have
density ~0.76/sqrt(ln x), and with K(n) ~ (log_4 n)^2/2 admissible pairs the
expected number of representations grows, so failures should be confined to
small n. The classical antecedent is Crocker's work on sums of two squares and
two powers of 2, where sparse summands likewise interact with the multiplicative
constraint of Fermat's theorem.

What that density argument misses is that the summands are not independent of
the constraint. Fermat's criterion is multiplicative and local, so an n whose
small-prime factorization is unfavourable — here n0 is divisible by 3, 7 and 19,
all == 3 (mod 4) — pushes many remainders simultaneously into residue classes
that cannot be sums of two squares. Lemma 4.1 makes this precise: 106 of the 136
remainders are excluded by congruences alone, before any arithmetic on n0 is
done. The counterexample is not a random coincidence among 136 independent
events but the product of one strong correlation.

### 6.2 Why certificates
The natural way to report a counterexample of this kind is "we ran an exhaustive
search and it found n0". That claim is only as good as the searcher — and, as
Section 5.6 shows, this project's own searcher had a real defect. Reducing the
result to 136 explicit prime powers decouples it from the search entirely: the
proof of Theorem 1 can be checked in 0.02 s by a program with no dependencies,
by a reader who trusts neither the searcher nor its author. We would encourage
the same separation for related conjectures in this family.

### 6.3 Beyond 2*10^9
Lemma 4.1 turns into a search strategy. Since a counterexample must have small
F(n), and F(n) is computable from n mod M by pure residue arithmetic without
touching a^2 + b^2, a search can filter candidates by F(n) first and run the
expensive coverage test only on the survivors. In the sampled range, requiring
F(n) <= 30 retains about 2*10^-5 of all n. The heuristic of Section 4.4 predicts
this will find nothing, which is precisely why it is worth running: it is a cheap
falsifiable prediction.

### 6.4 On AI assistance
The search, the certificate construction and this note were produced by the
author working with Anthropic's Claude, and this revision — including the
discovery of the v1.0.0 searcher defect — came out of the same workflow. We
state this plainly rather than in a footnote, because the more interesting point
is methodological: the reliability of the result does not depend on that
disclosure being taken on trust.

Theorem 1 reduces to a single integer, a finite list of 136 cases, and a
54-line standard-library checker. A reader who distrusts the provenance
entirely can verify it from those three objects alone, and the correction in
Section 5.6 is evidence the process is falsifiable in practice: a wrong comment
in the searcher survived the original release and was caught by re-derivation,
not by authority. Mathematical responsibility for the claims here rests with the
human author.
