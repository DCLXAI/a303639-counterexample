from math import comb
import json, subprocess, sys, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def test_pair_count():
    n = 800322180
    B = [comb(2*k+1,k) for k in range(16)]
    assert comb(33,16) > n
    assert sum(1 for i in range(16) for j in range(i,16) if B[i]+B[j] <= n) == 136

def test_certificates_valid():
    r = subprocess.run([sys.executable, os.path.join(ROOT,'certificates','check_certificates.py'),
                        os.path.join(ROOT,'certificates','certificates.json')],
                       capture_output=True, text=True)
    assert r.returncode == 0 and 'ALL 136' in r.stdout


def test_matches_official_oeis_data():
    """First 80 terms match the published DATA section of A303639."""
    from math import isqrt
    official=[0,1,1,2,1,3,2,2,1,2,3,3,3,3,4,2,2,2,3,4,4,5,2,4,1,2,3,3,5,3,5,1,3,1,1,6,3,8,3,6,
              2,4,4,2,7,5,6,2,5,2,4,5,4,8,4,7,2,4,1,3,6,4,7,3,5,2,4,2,4,9,5,6,2,6,4,5,4,7,5,2]
    B=[comb(2*k+1,k) for k in range(16)]
    def a(n):
        cnt=0
        for c in range(16):
            if B[c]>n: break
            for d in range(c,16):
                s=B[c]+B[d]
                if s>n: break
                r=n-s
                for x in range(isqrt(r//2)+1):
                    y2=r-x*x; y=isqrt(y2)
                    if y*y==y2 and y>=x: cnt+=1
        return cnt
    assert [a(n) for n in range(1,81)] == official
