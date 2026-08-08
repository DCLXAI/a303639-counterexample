#!/usr/bin/env python3
"""Why does the counterexample exist? Modular-obstruction analysis.

Naive heuristic: by Landau-Ramanujan the density of sums of two squares near
x is ~ 0.7642/sqrt(ln x) (~0.17 at x ~ 8*10^8). If the 136 remainders were
independent, P(all miss) ~ (1-0.17)^136 ~ e^-25, and one would expect ~0.03
counterexamples below 2*10^9 -- roughly consistent with finding exactly one,
but only because the misses are NOT independent:

n = 800322180 = 2^2 * 3 * 5 * 7 * 19 * 100291 is divisible by 3, 7 and 19
(all == 3 mod 4). This places the 136 remainders in correlated residue
classes: e.g. r == 3 or 6 (mod 9) forces v_3(r) = 1 (odd) outright. In the
event, small primes do most of the certifying work:
  p=3: 50/136,  p in {7,11,19,23}: 45/136,  large p: remaining 41.
Run this script to reproduce the table.
"""
from math import comb
from collections import Counter
import sympy

n = 800322180
B = [comb(2*k+1,k) for k in range(16)]
cert = Counter(); r9 = Counter()
for i in range(16):
    for j in range(i,16):
        s = B[i]+B[j]
        if s > n: break
        r = n-s
        bad = sorted(p for p,e in sympy.factorint(r).items() if p%4==3 and e%2==1)
        cert[bad[0]] += 1; r9[r%9] += 1
print("n =", n, "=", sympy.factorint(n))
print("smallest certifying prime distribution:", dict(sorted(cert.items())))
print("r mod 9 distribution:", dict(sorted(r9.items())),
      "\n(r == 3, 6 mod 9 => v_3(r) = 1, an automatic certificate)")
