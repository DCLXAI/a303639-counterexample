# A303639 Counterexample Project

**Claim.** `n = 800,322,180` is not of the form `a^2 + b^2 + C(2c+1,c) + C(2d+1,d)`
(`a<=b`, `c<=d`, all nonnegative), refuting the conjecture of Zhi-Wei Sun (2018)
recorded in [OEIS A303639](https://oeis.org/A303639) that `a(n) > 0` for all `n > 1`
(previously verified by the author for `n <= 6*10^8`). An exhaustive search shows
this is the **only** exception with `1 < n <= 2*10^9`.

## Structure
- `searcher/search.c` — self-contained searcher: 24.4 s / 500 MB single-core for the full `[1, 2*10^9]` scan (third independent implementation; re-finds the unique exception)
- `verifier/verify.py` — independent verification (sympy factorization; seconds)
- `certificates/certificates.json` — 136 machine-checkable certificates: for each
  qualifying pair `(c,d)` a prime `p == 3 (mod 4)` with odd exponent `e`, `p^e || r`
- `certificates/check_certificates.py` — trustless checker, **standard library only**
  (deterministic Miller-Rabin; no factorization needed) — completeness of the pair
  list is checked too
- `analysis/forced_table.py` — generates + verifies the Lemma 4.1 table: 106 of 136
  pairs are excluded by residues alone (`mod 8*9*49*361`)
- `analysis/obstruction.py`, `analysis/model.py` — heuristic model: obstruction
  amplification ~7*10^3; heuristically the unique counterexample, ever (NOTE.md §4)
- `tests/` — pytest sanity suite

## Reproduce
```
python3 certificates/check_certificates.py certificates/certificates.json
python3 verifier/verify.py          # requires sympy
python3 analysis/obstruction.py     # requires sympy
```

## Status
- [ ] OEIS comment submitted / accepted
- [ ] Author (Z.-W. Sun) notified
- [ ] Research note published

## Provenance
Search and analysis carried out by SU (Sunsu Jeong) with substantial assistance
from Anthropic's Claude; certificates independently machine-checked by the
standard-library verifier above.
