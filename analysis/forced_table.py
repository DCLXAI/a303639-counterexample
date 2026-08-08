#!/usr/bin/env python3
"""Generates and verifies the explicit table for Lemma 4.1 (deterministic,
no sampling): for n == 800322180 (mod M), M = 8*9*49*361 = 1273608, exactly
106 of the 136 qualifying pairs (c,d) have remainders that residues alone
prove are not sums of two squares. Output: ../certificates/forced_table.json
"""
import json, os
from math import comb

n0, M = 800322180, 8*9*49*361
B = [comb(2*k+1,k) for k in range(16)]

def classify(rho):
    """rho = r mod M. Returns obstruction case or None. Each case is a
    rigorous implication for EVERY r == rho (mod M):
      mod8 : odd part of r == 3 (mod 4)  -> some p==3(4) to odd order
      v3   : r == 3,6 (mod 9)   -> v_3(r) = 1
      v7   : 7|r, 49!|r (mod 49)  -> v_7(r) = 1
      v19  : 19|r, 361!|r (mod 361) -> v_19(r) = 1
    """
    r8, r9, r49, r361 = rho%8, rho%9, rho%49, rho%361
    if r8 in (3,7) or r8 == 6:                 return 'mod8'
    if r9 in (3,6):                            return 'v3'
    if r49 % 7 == 0 and r49 != 0:              return 'v7'
    if r361 % 19 == 0 and r361 != 0:           return 'v19'
    return None

rows, free = [], []
for i in range(16):
    for j in range(i,16):
        if B[i]+B[j] > n0: continue
        rho = (n0 - B[i] - B[j]) % M
        case = classify(rho)
        (rows if case else free).append(dict(c=i, d=j, r_mod_M=rho, case=case))

# sanity: counts must match Lemma 4.1
from collections import Counter
cnt = Counter(r['case'] for r in rows)
assert len(rows) == 106 and len(free) == 30, (len(rows), len(free))
assert cnt == {'mod8':51, 'v3':29, 'v7':13, 'v19':13}, cnt
# soundness spot-check on actual n0 remainders
for r in rows:
    rr = n0 - B[r['c']] - B[r['d']]
    assert rr % M == r['r_mod_M']
    if r['case']=='v3':  assert rr%3==0 and rr%9!=0
    if r['case']=='v7':  assert rr%7==0 and rr%49!=0
    if r['case']=='v19': assert rr%19==0 and rr%361!=0
    if r['case']=='mod8':
        t=rr
        while t%2==0: t//=2
        assert t%4==3

out = os.path.join(os.path.dirname(__file__), '..', 'certificates', 'forced_table.json')
json.dump(dict(n0=n0, M=M, forced=rows, free=[dict(c=f['c'],d=f['d']) for f in free],
               counts=dict(cnt)), open(out,'w'), indent=1)
print(f"Lemma 4.1 table verified and written: 106 forced ({dict(cnt)}), 30 free.")
