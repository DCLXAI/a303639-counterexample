/* A303639 counterexample searcher.
 * Finds all n in (1, N] with no representation n = a^2 + b^2 + B(c) + B(d),
 * where B(k) = binomial(2k+1, k), a<=b, c<=d, all nonnegative.
 * Method: bitset coverage. base = { a^2 + b^2 <= N }, then OR of 136 shifted
 * copies (shifts B(c)+B(d)); uncovered n are counterexamples.
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
    /* B(k) = C(2k+1,k), k=0..15 (B(15) <= 2e9 < B(16)) */
    long long B[16]; B[0]=1;
    for(int k=1;k<16;k++){                 /* B(k) = B(k-1)*(2k)*(2k+1)/(k*(k+1)) */
        __int128 t=(__int128)B[k-1]*(2*k)*(2*k+1)/((long long)k*(k+1));
        B[k]=(long long)t;
    }
    /* 136 shifts: cov |= base << (B[i]+B[j]) */
    int shifts=0;
    for(int i=0;i<16;i++) for(int j=i;j<16;j++){
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
