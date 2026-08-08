#!/usr/bin/env python3
"""Machine-checkable certificate verifier for the A303639 counterexample.

Checks, WITHOUT any factorization and using only the standard library, that
for every pair (c,d) with B(c)+B(d) <= n  (B(k) = binomial(2k+1,k)):
  1. the pair list is complete (all 136 qualifying pairs present, none missed),
  2. r = n - B(c) - B(d) is correct,
  3. the certificate prime p is prime (deterministic Miller-Rabin, 64-bit),
  4. p == 3 (mod 4),  e is odd,  and p^e || r  (p^e divides r, p^(e+1) does not).
By Fermat's two-squares theorem each certified r is not a sum of two squares,
hence n has no representation a^2+b^2+B(c)+B(d) and A303639(n) = 0.
"""
import json, sys
from math import comb

def is_prime(m):
    if m < 2: return False
    for q in (2,3,5,7,11,13,17,19,23,29,31,37):
        if m % q == 0: return m == q
    d, s = m-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):   # 결정적 (m < 3.3*10^24)
        x = pow(a, d, m)
        if x in (1, m-1): continue
        for _ in range(s-1):
            x = x*x % m
            if x == m-1: break
        else: return False
    return True

data = json.load(open(sys.argv[1] if len(sys.argv)>1 else 'certificates.json'))
n = data['n']; certs = {(c['c'],c['d']): c for c in data['certificates']}
B = [comb(2*k+1,k) for k in range(16)]
assert comb(33,16) > n, "B(16) bound"
expected = [(i,j) for i in range(16) for j in range(i,16) if B[i]+B[j] <= n]
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
