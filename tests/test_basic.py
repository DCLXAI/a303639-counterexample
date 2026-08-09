from math import comb
import json, shutil, subprocess, sys, os, tempfile
import pytest
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def kmax(N, B):
    return max(k for k in range(len(B)) if B[k] <= N)

def npairs(N, B):
    K = kmax(N, B)
    return sum(1 for i in range(K+1) for j in range(i,K+1) if B[i]+B[j] <= N)

def test_pair_count():
    """At n0 the cutoff really is k <= 15, because B(16) exceeds n0."""
    n = 800322180
    B = [comb(2*k+1,k) for k in range(20)]
    assert comb(33,16) > n
    assert kmax(n, B) == 15
    assert npairs(n, B) == 136

def test_index_16_qualifies_below_2e9():
    """Regression: B(16) < 2e9, so a full [1,2e9] scan needs 152 pairs, not 136.

    Releases up to v1.0.0 hard-coded k <= 15 in the searcher on the strength of
    the false claim "B(15) <= 2e9 < B(16)". B(16) is 1166803110 -- comfortably
    below 2e9 -- so 16 pairs were missing from the scan.
    """
    B = [comb(2*k+1,k) for k in range(20)]
    assert B[16] == 1166803110 < 2*10**9 < B[17]
    assert kmax(2*10**9, B) == 16
    assert npairs(2*10**9, B) == 152
    # first n at which an index-16 pair becomes available
    assert B[16] + B[0] == 1166803111
    assert npairs(B[16]+B[0], B) == 137 and npairs(B[16]+B[0]-1, B) == 136

def test_certificates_valid():
    r = subprocess.run([sys.executable, os.path.join(ROOT,'certificates','check_certificates.py'),
                        os.path.join(ROOT,'certificates','certificates.json')],
                       capture_output=True, text=True)
    assert r.returncode == 0 and 'ALL 136' in r.stdout


def a303639(n):
    """Direct definition of a(n), used to cross-check against the OEIS entry."""
    from math import isqrt
    B = [comb(2*k+1,k) for k in range(20)]
    cnt = 0
    for c in range(20):
        if B[c] > n: break
        for d in range(c,20):
            s = B[c]+B[d]
            if s > n: break
            r = n-s
            for x in range(isqrt(r//2)+1):
                y2 = r-x*x; y = isqrt(y2)
                if y*y == y2 and y >= x: cnt += 1
    return cnt


def test_matches_official_oeis_data():
    """First 80 terms match the published DATA section of A303639."""
    official=[0,1,1,2,1,3,2,2,1,2,3,3,3,3,4,2,2,2,3,4,4,5,2,4,1,2,3,3,5,3,5,1,3,1,1,6,3,8,3,6,
              2,4,4,2,7,5,6,2,5,2,4,5,4,8,4,7,2,4,1,3,6,4,7,3,5,2,4,2,4,9,5,6,2,6,4,5,4,7,5,2]
    assert [a303639(n) for n in range(1,81)] == official


def test_matches_oeis_example_values():
    """The n with a(n) = 1 listed in the entry's EXAMPLE section."""
    assert all(a303639(n) == 1 for n in (2530, 3258, 5300, 13453, 20964))


@pytest.mark.skipif(not (shutil.which('gcc') or shutil.which('cc')),
                    reason='no C compiler available')
def test_searcher_includes_index_16():
    """The searcher must derive its index cutoff from N, not hard-code it.

    Run at N = B(16)+B(0) = 1166803111, the smallest bound at which an
    index-16 pair qualifies: a correct searcher reports k = 0..16 and 137
    shifts. The pre-v1.0.1 searcher reports 136 and never touches k = 16.
    Costs ~300 MB and a few seconds; the old N=1e8 smoke test could not
    detect this, since B(15) alone already exceeds 1e8.
    """
    cc = shutil.which('gcc') or shutil.which('cc')
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, 'search')
        build = subprocess.run([cc, '-O3', '-o', exe, os.path.join(ROOT,'searcher','search.c')],
                               capture_output=True, text=True)
        assert build.returncode == 0, build.stderr
        r = subprocess.run([exe, '1166803111'], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        assert 'k = 0..16' in r.stderr, r.stderr
        assert '137 shifts done' in r.stderr, r.stderr
        # n0 < B(16), so it stays the only uncovered value at this bound too
        assert r.stdout.strip() == 'COUNTEREXAMPLE n=800322180'
