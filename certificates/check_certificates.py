#!/usr/bin/env python3
"""Machine-checkable certificate verifier for the A303639 counterexample.

Checks, WITHOUT any factorization and using only the standard library, that
for every pair (c,d) with B(c)+B(d) <= n  (B(k) = binomial(2k+1,k)):
  1. the pair list is complete (all 136 qualifying pairs present, none missed),
  2. r = n - B(c) - B(d) is correct,
  3. the certificate prime p is prime (deterministic trial division),
  4. p == 3 (mod 4),  e is odd,  and p^e || r  (p^e divides r, p^(e+1) does not).
By Fermat's two-squares theorem each certified r is not a sum of two squares,
hence n has no representation a^2+b^2+B(c)+B(d) and A303639(n) = 0.
"""
import json, sys
from math import comb

def is_prime(m):
    """Unconditional trial division -- deliberately the dumbest possible test.

    Every certificate prime here satisfies p <= r < n < 2^31, so isqrt(p) <
    46341 and the whole check costs a few hundred thousand operations. Using
    Miller-Rabin instead would make correctness depend on a published bound on
    the smallest strong pseudoprime for the chosen base set, which a referee
    would then have to look up and audit; trial division needs no such input.
    """
    if m < 2: return False
    if m < 4: return True
    if m % 2 == 0: return False
    f = 3
    while f*f <= m:
        if m % f == 0: return False
        f += 2
    return True

data = json.load(open(sys.argv[1] if len(sys.argv)>1 else 'certificates.json'))
n = data['n']; certs = {(c['c'],c['d']): c for c in data['certificates']}
# Removes the "k <= 15" assumption: derive the index cutoff from n itself.
# B(k) is increasing, so the largest usable index is the last k with B(k) <= n.
B, k = [1], 0
while B[k] <= n:
    k += 1; B.append(comb(2*k+1, k))
KMAX = k-1                                     # B(KMAX) <= n < B(KMAX+1)
assert B[KMAX] <= n < B[KMAX+1]
assert KMAX == 15 and B[16] == 1166803110, (KMAX, B[16])   # for n = 800322180
expected = [(i,j) for i in range(KMAX+1) for j in range(i,KMAX+1) if B[i]+B[j] <= n]
assert set(certs) == set(expected), "pair list incomplete or has extras"
for (i,j) in expected:
    c = certs[(i,j)]
    assert c['Bc']==B[i] and c['Bd']==B[j], (i,j)
    r = n - B[i] - B[j]
    assert c['r']==r and r>0, (i,j)
    p, e = c['p'], c['e']
    assert p % 4 == 3 and e % 2 == 1 and is_prime(p), (i,j,p,e)
    assert r % p**e == 0 and r % p**(e+1) != 0, (i,j,p,e)
print(f"ALL {len(expected)} CERTIFICATES VALID -> A303639({n}) = 0. Conjecture refuted.")
