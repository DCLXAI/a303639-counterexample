#!/usr/bin/env python3
"""Quantitative obstruction model for the A303639 counterexample.
Reproduces (seeded) the numbers cited in NOTE.md Section 4. Stdlib + numpy.
Runtime: ~2-3 min."""
import numpy as np, random
from math import comb, log, sqrt, exp

n0   = 800322180
M    = 8*9*49*361                    # obstruction modulus
NMAX = 2*10**9                       # upper end of the sampling range

# The pair list must be derived from the sampling range, not from n0. B(16) =
# 1166803110 < NMAX, so every n >= B(16)+B(0) = 1166803111 admits pairs with
# d = 16 that do not qualify at n0. Fixing the list at n0's 136 pairs would
# silently measure a restricted problem for the upper half of the range.
B, k = [1], 0
while B[k] <= NMAX:
    k += 1; B.append(comb(2*k+1, k))
KMAX = k-1                                                  # = 16 at NMAX = 2e9
ALL = [(i,j) for i in range(KMAX+1) for j in range(i,KMAX+1) if B[i]+B[j] <= NMAX]
assert KMAX == 16 and len(ALL) == 152, (KMAX, len(ALL))

def qual(n):
    """Pairs (c,d), c<=d, actually available at n."""
    return [(i,j) for i,j in ALL if B[i]+B[j] <= n]

pairs = qual(n0)                                            # 136 for n0
assert len(pairs) == 136

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

# 2. free-pair count F(n) and total qualifying-pair count K(n) for random n
#    (pure modular arithmetic).  K(n) is recorded too: it jumps from 136 to 152
#    across the range as the index-16 pairs come into play.
Fs=[]; Ks=[]
for _ in range(300000):
    n=random.randrange(4*10**8, NMAX)
    q=[(i,j) for i,j in ALL if B[i]+B[j]<=n]
    Fs.append(sum(1 for i,j in q if forced_fail(n-B[i]-B[j]) is None)); Ks.append(len(q))
Fs=np.array(Fs); Ks=np.array(Ks)
print(f"[2] F(n): mean {Fs.mean():.2f}, min {Fs.min()};  n0 has F={free} "
      f"-> percentile {100*(Fs<=free).mean():.5f}%")
print(f"    K(n): mean {Ks.mean():.1f}, range {Ks.min()}..{Ks.max()}")

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
    qp=qual(n)
    Fn=sum(1 for i,j in qp if forced_fail(n-B[i]-B[j]) is None)
    obs.append((sum(not is_S2(n-B[i]-B[j]) for i,j in qp), len(qp), Fn))
fr=np.array([a for a,_,_ in obs]); tt=np.array([b for _,b,_ in obs])
fn=np.array([c for _,_,c in obs])
pred=(tt-fn)+(1-d_cond)*fn          # forced pairs fail w.p. 1, free ones w.p. 1-d_cond
print(f"[4] same-class n (300): observed mean fails {fr.mean():.1f}/{tt.mean():.1f}, "
      f"max {fr.max()}; full exceptions {(fr==tt).sum()}")
print(f"    model predicts mean fails {pred.mean():.1f} (observed {fr.mean():.1f})")

# 5. expectations
d_marg=0.174                        # measured marginal density of sums of two squares
P_naive=(1-d_marg)**136             # n0 has 136 qualifying pairs
P_cond=(1-d_cond)**free             # ... of which only `free` are not forced
print(f"[5] n0: P(all 136 fail) naive ~ {P_naive:.2e}; conditioned on its class "
      f"~ {P_cond:.2e}  -> amplification ~ {P_cond/P_naive:.3g}")
E_naive=np.mean((1-d_marg)**Ks.astype(float))*(2*10**9-4*10**8)
E=np.mean((1-d_cond)**Fs.astype(float))*(2*10**9-4*10**8)
print(f"[5] E[#exceptions, 4e8..2e9]: naive ~ {E_naive:.3f}, refined ~ {E:.3f}")
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
