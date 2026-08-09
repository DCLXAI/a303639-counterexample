/* A303639 counterexample searcher.
 * Finds all n in (1, N] with no representation n = a^2 + b^2 + B(c) + B(d),
 * where B(k) = binomial(2k+1, k), a<=b, c<=d, all nonnegative.
 * Method: bitset coverage. base = { a^2 + b^2 <= N }, then OR of one shifted
 * copy per qualifying pair (shifts B(c)+B(d)); uncovered n are counterexamples.
 * At N = 2e9 there are 152 such pairs (B(16) = 1166803110 < 2e9).
 * Build: gcc -O3 -march=native -o search search.c
 * Run:   ./search 2000000000        (~250 MB x2 memory)
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define SETB(b,v) (b[(v)>>6]|=1ULL<<((v)&63))
#define GETB(b,v) ((b[(v)>>6]>>((v)&63))&1ULL)
static double el(clock_t t0){ return (double)(clock()-t0)/CLOCKS_PER_SEC; }
int main(int argc, char **argv){
    long long N = argc>1 ? atoll(argv[1]) : 2000000000LL;
    long long W = N/64 + 2;
    clock_t t0 = clock();
    uint64_t *base = calloc(W,8), *cov = calloc(W,8);
    if(!base || !cov){ fprintf(stderr,"alloc fail\n"); return 1; }
    /* base: sums of two squares (direct double loop, ~pi/8*N iterations) */
    for(long long a=0; a*a<=N; a++){
        long long aa=a*a;
        for(long long b=a; aa+b*b<=N; b++) SETB(base, aa+b*b);
    }
    fprintf(stderr,"base fill: %.1fs\n", el(t0));
    /* B(k) = C(2k+1,k), generated while B(k) <= N -- do not hard-code the
     * cutoff: B(16) = 1166803110 is itself below 2e9, so k = 16 qualifies at
     * the full search bound. B grows like 4^k, so K stays under 32 for any
     * N a long long can hold. */
    long long B[64]; B[0]=1; int K=0;
    for(;;){                               /* B(k) = B(k-1)*(2k)*(2k+1)/(k*(k+1)) */
        int k=K+1;
        __int128 t=(__int128)B[K]*(2*k)*(2*k+1)/((long long)k*(k+1));
        if(t>(__int128)N) break;
        B[++K]=(long long)t;
    }
    fprintf(stderr,"B(k) <= N for k = 0..%d (B(%d) = %lld)\n", K, K, B[K]);
    /* one shift per qualifying pair: cov |= base << (B[i]+B[j]) */
    int shifts=0;
    for(int i=0;i<=K;i++) for(int j=i;j<=K;j++){
        long long s=B[i]+B[j];
        if(s>N) continue;
        long long wq=s>>6; int r=s&63;
        if(r==0){ for(long long w=W-1; w>=wq; w--) cov[w]|=base[w-wq]; }
        else { for(long long w=W-1; w>wq; w--)
                   cov[w] |= (base[w-wq]<<r) | (base[w-wq-1]>>(64-r));
               cov[wq]|=base[0]<<r; }
        shifts++;
    }
    fprintf(stderr,"%d shifts done: %.1fs\n", shifts, el(t0));
    long long found=0;
    for(long long n=2;n<=N;n++)
        if(!GETB(cov,n)){ printf("COUNTEREXAMPLE n=%lld\n", n); found++; }
    fprintf(stderr,"scan done: %.1fs total, %lld counterexample(s), mem %lld MB\n",
            el(t0), found, 2*W*8/1000000);
    return 0;
}
