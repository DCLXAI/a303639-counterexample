# A303639 Counterexample Project

[![DOI](https://zenodo.org/badge/1328043397.svg)](https://doi.org/10.5281/zenodo.21863025)

**Claim.** `n = 800,322,180` is not of the form `a^2 + b^2 + C(2c+1,c) + C(2d+1,d)`
(`a<=b`, `c<=d`, all nonnegative), refuting the conjecture of Zhi-Wei Sun (2018)
recorded in [OEIS A303639](https://oeis.org/A303639) that `a(n) > 0` for all `n > 1`
(previously verified by the author for `n <= 6*10^8`). An exhaustive search shows
this is the **only** exception with `1 < n <= 2*10^9`.

The proof of the claim is finite and self-contained: 136 explicit prime-power
certificates, checkable in 0.02 s by a standard-library-only program that trusts
neither the searcher nor its author.

## Structure
- `searcher/search.c` — self-contained searcher: 6.1 s / 500 MB single-core for the
  full `[1, 2*10^9]` scan on an Apple M4 (third independent implementation;
  re-finds the unique exception). The scan uses **152** shifts: `B(16) = 1166803110`
  is below `2*10^9`, even though it exceeds `n` itself
- `verifier/verify.py` — independent verification (sympy factorization; seconds)
- `certificates/certificates.json` — 136 machine-checkable certificates: for each
  qualifying pair `(c,d)` a prime `p == 3 (mod 4)` with odd exponent `e`, `p^e || r`
- `certificates/check_certificates.py` — trustless checker, **standard library only**
  (unconditional trial division; no factorization and no primality-test bound to
  audit) — the pair list is re-derived from `n`, so completeness is checked too
- `analysis/forced_table.py` — generates + verifies the Lemma 4.1 table: 106 of 136
  pairs are excluded by residues alone (`mod 8*9*49*361`)
- `analysis/obstruction.py`, `analysis/model.py` — heuristic model: obstruction
  amplification ~1.3*10^4; heuristically the unique counterexample, ever (NOTE.md §4)
- `tests/` — pytest suite, including a regression test for the `B(16)` cutoff

## Reproduce
```
python3 certificates/check_certificates.py certificates/certificates.json
python3 verifier/verify.py          # requires sympy
python3 analysis/forced_table.py
python3 analysis/model.py           # requires numpy + sympy; ~2-3 min
pytest tests/ -q
gcc -O3 -o search searcher/search.c && ./search 2000000000
```

## Changes since v1.0.0
v1.0.0 contained four defects, none affecting the counterexample itself
(its certificates were valid as published) but two affecting the code:
the searcher hard-coded a pair cutoff of `k <= 15` on a false premise and so ran
136 shifts instead of 152; `analysis/model.py` fixed its pair list at `n`'s 136
pairs while sampling up to `2*10^9`, so the §4.2-4.4 statistics measured a
restricted problem; Lemma 4.1 was stated without the range bound it needs; and
the checker's Miller-Rabin bases carried an incorrect deterministic bound.
All are corrected and the §4 figures recomputed — see **NOTE.md §5.6**.
Cite the corrected release — v1.0.1, DOI [10.5281/zenodo.21863408](https://doi.org/10.5281/zenodo.21863408)
— not v1.0.0 (DOI 10.5281/zenodo.21863026). The badge above is the *concept*
DOI, which always resolves to the newest version.

## Status
- [ ] OEIS comment submitted / accepted
- [ ] Author (Z.-W. Sun) notified
- [ ] Research note published

## Provenance
Search and analysis carried out by SU (Sunsu Jeong) with substantial assistance
from Anthropic's Claude; certificates independently machine-checked by the
standard-library verifier above. Mathematical responsibility rests with the
human author (NOTE.md §6.4).
