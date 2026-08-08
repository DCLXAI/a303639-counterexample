#!/usr/bin/env python3
"""
Independent verification: n = 800322180 is NOT of the form
a^2 + b^2 + binomial(2c+1,c) + binomial(2d+1,d)  (a<=b, c<=d, all >= 0),
refuting Zhi-Wei Sun's conjecture in OEIS A303639 (a(n) > 0 for all n > 1,
previously verified for n = 2..6*10^8).

Method: for every pair (c,d) with B(c)+B(d) <= n, the remainder
r = n - B(c) - B(d) must fail the sum-of-two-squares criterion,
i.e., r has a prime p == 3 (mod 4) with odd exponent (Fermat).
Requires: sympy  (pip install sympy)
"""
from math import comb
import sympy

n = 800322180
B = [comb(2*k+1, k) for k in range(16)]      # B(15)=300540195 <= n < B(16)
assert comb(33, 16) > n

pairs = 0
for i in range(16):
    for j in range(i, 16):
        s = B[i] + B[j]
        if s > n:
            break
        r = n - s
        f = sympy.factorint(r)
        bad = [(p, e) for p, e in f.items() if p % 4 == 3 and e % 2 == 1]
        assert bad, f"FAIL: r={r} (c={i}, d={j}) IS a sum of two squares!"
        pairs += 1

print(f"VERIFIED: all {pairs} pairs (c,d) leave a non-two-square remainder.")
print("Therefore a(800322180) = 0 and the A303639 conjecture is false.")
