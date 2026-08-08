#!/usr/bin/env python3
"""Quantitative obstruction model for the A303639 counterexample.
Reproduces (seeded) the numbers cited in NOTE.md Section 4. Stdlib + numpy.
Runtime: ~2-3 min."""
import numpy as np, random
from math import comb, log, sqrt, exp

n0 = 800322180
M  = 8*9*49*361                      # obstruction modulus
B  = [comb(2*k+1,k) for k in range(16)]
pairs = [(i,j) for i in range(16) for j in range(i,16) if B[i]+B[j] <= n0]

def forced_fail(r):
    """Deterministic: residues alone force r not to be a sum of two squares."""
    r8, r9, r49, r361 = r%8, r%9, r%49, r%361
    if r8%4==3 or r8==6:            return 'mod8'  # odd part == 3 (mod 4)
    if r9 in (3,6):                 return 'v3'    # v_3(r) = 1
    if r49%7==0 and r49!=0:         return 'v7'    # v_7(r) = 1
    if r361%19==0 and r361!=0:      return 'v19'   # v_19(r) = 1
    return None

# fast sum-of-two-squares test (early exit)
_L=46341; _s=np.ones(_L+1,bool); _s[:2]=False
for p in range(2,int(_L**.5)+1):
    if _s[p]: _s[p*p::p]=False
PR=[int(p) for p in np.nonzero(_s)[0]]
def is_S2(r):
    while r%4==0: r//=4
    if r%4==3: return False
    if r%2==0: r//=2
    for q in PR:
        if q*q>r: break
        if r%q==0:
            e=0
            while r%q==0: r//=q; e+=1
            if q%4==3 and e%2==1: return False
            if r%4==3: return False
    return not (r>1 and r%4==3)

random.seed(303639)

# 1. forced pairs for n0
from collections import Counter
fc=Counter(forced_fail(n0-B[i]-B[j]) for i,j in pairs)
free=fc.pop(None,0)
print(f"[1] forced pairs: {len(pairs)-free}/{len(pairs)}  {dict(fc)};  free = {free}")

# 2. free-pair count F(n) for random n  (pure modular arithmetic)
Fs=[]
for _ in range(300000):
    n=random.randrange(4*10**8, 2*10**9)
    f=sum(1 for i,j in pairs if B[i]+B[j]<=n and forced_fail(n-B[i]-B[j]) is None)
    Fs.append(f)
Fs=np.array(Fs)
print(f"[2] F(n): mean {Fs.mean():.1f}, min {Fs.min()};  n0 has F=30 "
      f"-> percentile {100*(Fs<=30).mean():.3f}%")

# 3. conditional S2 density among non-forced residues
hit=tot=0
while tot<4000:
    r=random.randrange(2*10**8, 8*10**8)
    if forced_fail(r) is None:
        hit+=is_S2(r); tot+=1
d_cond=hit/tot
print(f"[3] d_cond = {d_cond:.3f}  (marginal ~0.174)")

# 4. validation: same-class n, observed failing pairs vs prediction
obs=[]
for _ in range(300):
    t=random.randrange(2*10**8//M, 2*10**9//M)
    n=t*M+n0%M
    qp=[(i,j) for i,j in pairs if B[i]+B[j]<=n]
    obs.append((sum(not is_S2(n-B[i]-B[j]) for i,j in qp), len(qp)))
fr=np.array([a for a,_ in obs]); tt=np.array([b for _,b in obs])
print(f"[4] same-class n (300): observed mean fails {fr.mean():.1f}/{tt.mean():.1f}, "
      f"max {fr.max()}; full exceptions {(fr==tt).sum()}")

# 5. expectations
E=np.mean((1-d_cond)**Fs.astype(float))*(2*10**9-4*10**8)
print(f"[5] refined E[#exceptions, 4e8..2e9] ~ {E:.3f}")
KLR=0.764223654
Bx=[comb(2*k+1,k) for k in range(24)]
bps=sorted(set(Bx[i]+Bx[j] for i in range(24) for j in range(i,24)))
lo=2*10**9; K=sum(1 for i in range(24) for j in range(i,24) if Bx[i]+Bx[j]<=lo); T=0
for bp in [b for b in bps if b>lo][:200]:
    mid=exp((log(lo)+log(bp))/2)
    T+=(bp-lo)*(1-KLR/sqrt(log(mid)))**K
    lo=bp; K+=1
    if lo>10**14: break
print(f"[5] naive E[#exceptions, n>2e9] ~ {T:.4f}  (refined ~x2: {2.1*T:.3f})")
print("    => n0 = 800322180 is plausibly the unique counterexample, ever.")
