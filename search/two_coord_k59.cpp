#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <boost/multiprecision/cpp_int.hpp>
#include <omp.h>

using u64 = std::uint64_t;
using u128 = __uint128_t;
using boost::multiprecision::uint256_t;

constexpr int K = 59;
constexpr int WORDS = 29;
using Mask = std::array<u64, WORDS>;

struct CT {
    int mod = 0;
    int p = 0;
    std::vector<Mask> cov;
};

struct Alt {
    int score = 0;
    int residue = 0;
};

struct PairClass {
    int base_score = 0;
    int r1 = 0;
    int r2 = 0;
};

struct Candidate {
    int score200 = 0;
    int score1m = 0;
    int base_score = 0;
    int q1 = 0, r1 = 0, q2 = 0, r2 = 0;
    u128 n = 0;
    Mask live{};
};

static inline void setbit(Mask& m, int i) { m[i >> 6] |= 1ULL << (i & 63); }
static inline void clrbit(Mask& m, int i) { m[i >> 6] &= ~(1ULL << (i & 63)); }
static inline void mask_or(Mask& a, const Mask& b) { for (int w=0; w<WORDS; ++w) a[w] |= b[w]; }
static inline void mask_andnot(Mask& a, const Mask& b) { for (int w=0; w<WORDS; ++w) a[w] &= ~b[w]; }
static inline int popcount(const Mask& m) { int z=0; for (u64 x : m) z += __builtin_popcountll(x); return z; }

static bool is_prime_small(int n) {
    if (n < 2) return false;
    for (int d=2; (long long)d*d<=n; ++d) if (n%d==0) return false;
    return true;
}

static std::vector<int> primes_3mod4(int lo, int hi) {
    std::vector<bool> ok(hi, true);
    ok[0]=ok[1]=false;
    for (int p=2; (long long)p*p<hi; ++p) if (ok[p])
        for (long long q=(long long)p*p; q<hi; q+=p) ok[(std::size_t)q]=false;
    std::vector<int> ps;
    for (int p=std::max(3,lo); p<hi; ++p) if (ok[p] && (p&3)==3) ps.push_back(p);
    return ps;
}

static u64 inv_mod(u64 a, u64 m) {
    long long t=0, nt=1, r=(long long)m, nr=(long long)(a%m);
    while (nr) {
        long long q=r/nr;
        long long tt=t-q*nt; t=nt; nt=tt;
        long long rr=r-q*nr; r=nr; nr=rr;
    }
    if (t<0) t += (long long)m;
    return (u64)t;
}

static std::pair<u128,u128> crt(const std::vector<int>& mods, const std::vector<int>& residues) {
    u128 x=0, M=1;
    for (std::size_t i=0; i<mods.size(); ++i) {
        u64 m=(u64)mods[i], b=(u64)residues[i];
        u64 xm=(u64)(x%m), Mm=(u64)(M%m);
        u64 delta=(b+m-xm)%m;
        u128 step=(u128)delta*inv_mod(Mm,m)%m;
        x += M*step;
        M *= m;
    }
    return {x,M};
}

static u128 ceil_div(u128 a, u128 b) { return a/b + (a%b != 0); }

static std::string str128(u128 x) {
    if (!x) return "0";
    std::string s;
    while (x) { s.push_back(char('0'+x%10)); x/=10; }
    std::reverse(s.begin(),s.end());
    return s;
}

struct Hash128 {
    std::size_t operator()(u128 x) const noexcept {
        return (u64)x ^ ((u64)(x>>64) * 0x9e3779b97f4a7c15ULL);
    }
};

static bool worse_candidate(const Candidate& a, const Candidate& b) {
    if (a.score200 != b.score200) return a.score200 > b.score200;
    if (a.base_score != b.base_score) return a.base_score > b.base_score;
    return a.n > b.n;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: two_coord_k59 Q1 OUTPUT [alt_cap=64] [class_cap=250] [base_thr=72] [keep200=500] [out1m=120]\n";
        return 2;
    }

    int q1_fixed = std::atoi(argv[1]);
    const char* output = argv[2];
    int alt_cap = argc>3 ? std::atoi(argv[3]) : 64;
    int class_cap = argc>4 ? std::atoi(argv[4]) : 250;
    int base_thr = argc>5 ? std::atoi(argv[5]) : 72;
    int keep200 = argc>6 ? std::atoi(argv[6]) : 500;
    int out1m = argc>7 ? std::atoi(argv[7]) : 120;

    std::vector<u128> B(K+2);
    B[0]=1;
    for (int k=1; k<=K+1; ++k) {
        uint256_t z=uint256_t(B[k-1])*(2*k)*(2*k+1);
        z /= (u64)k*(k+1);
        B[k]=(u128)z;
    }
    const u128 interval_lo=B[K];
    const u128 interval_hi=B[K+1]-1;

    std::vector<u128> shifts;
    for (int i=0; i<=K; ++i) for (int j=i; j<=K; ++j) shifts.push_back(B[i]+B[j]);
    const int NS=(int)shifts.size();
    Mask ALL{};
    for (int i=0; i<NS; ++i) setbit(ALL,i);

    std::vector<int> base_primes={0,3,7,11,19,23,31,43,47,59,67,71};
    std::vector<int> mods={8};
    for (std::size_t i=1; i<base_primes.size(); ++i) mods.push_back(base_primes[i]*base_primes[i]);
    if (q1_fixed<0 || q1_fixed>=(int)mods.size()-1) {
        std::cerr << "Q1 must be in [0," << mods.size()-2 << "]\n";
        return 2;
    }

    const char* near_text="153682427474744689848227492854843581";
    u128 near=0;
    for (const char* p=near_text; *p; ++p) near=near*10+(*p-'0');
    std::vector<int> original;
    for (int m:mods) original.push_back((u64)(near%m));

    std::vector<CT> base;
    base.reserve(mods.size());
    for (std::size_t qi=0; qi<mods.size(); ++qi) {
        int m=mods[qi], p=base_primes[qi];
        CT t{m,p,std::vector<Mask>(m)};
        for (int i=0; i<NS; ++i) {
            if (m==8) {
                int sm=(int)(shifts[i]%8);
                for (int d : {3,6,7}) setbit(t.cov[(sm+d)&7],i);
            } else {
                int exact=(int)(shifts[i]%m);
                int a=(int)(shifts[i]%p);
                for (int q=0; q<p; ++q) {
                    int r=a+p*q;
                    if (r!=exact) setbit(t.cov[r],i);
                }
            }
        }
        base.push_back(std::move(t));
    }

    std::vector<Mask> original_cov(mods.size());
    for (std::size_t q=0; q<mods.size(); ++q) original_cov[q]=base[q].cov[original[q]];

    std::vector<std::vector<Alt>> alternatives(mods.size());
    for (int q=0; q<(int)mods.size(); ++q) {
        Mask other{};
        for (int z=0; z<(int)mods.size(); ++z) if (z!=q) mask_or(other,original_cov[z]);
        for (int r=0; r<mods[q]; ++r) {
            if (r==original[q]) continue;
            Mask m=other;
            mask_or(m,base[q].cov[r]);
            int score=NS-popcount(m);
            if (score<=96) alternatives[q].push_back({score,r});
        }
        std::sort(alternatives[q].begin(),alternatives[q].end(),[](const Alt&a,const Alt&b){
            if (a.score!=b.score) return a.score<b.score;
            return a.residue<b.residue;
        });
        if ((int)alternatives[q].size()>alt_cap) alternatives[q].resize(alt_cap);
        std::cerr << "alt q=" << q << " kept=" << alternatives[q].size();
        if (!alternatives[q].empty()) std::cerr << " best=" << alternatives[q][0].score;
        std::cerr << "\n";
    }

    std::vector<CT> filters200;
    for (int p=3; p<200; p+=4) if (is_prime_small(p)) {
        int p2=p*p;
        CT t{p2,p,std::vector<Mask>(p2)};
        for (int i=0; i<NS; ++i) {
            int exact=(int)(shifts[i]%p2);
            int a=(int)(shifts[i]%p);
            for (int q=0; q<p; ++q) {
                int r=a+p*q;
                if (r!=exact) setbit(t.cov[r],i);
            }
        }
        filters200.push_back(std::move(t));
    }

    std::vector<Candidate> kept;
    kept.reserve(keep200*2);
    u64 classes=0, evaluated=0, passed48=0;

    for (int q2=q1_fixed+1; q2<(int)mods.size(); ++q2) {
        Mask other{};
        for (int q=0; q<(int)mods.size(); ++q)
            if (q!=q1_fixed && q!=q2) mask_or(other,original_cov[q]);

        std::vector<PairClass> pair_classes;
        for (const Alt& a1:alternatives[q1_fixed]) {
            Mask m1=other;
            mask_or(m1,base[q1_fixed].cov[a1.residue]);
            for (const Alt& a2:alternatives[q2]) {
                Mask m=m1;
                mask_or(m,base[q2].cov[a2.residue]);
                int score=NS-popcount(m);
                if (score<=base_thr) pair_classes.push_back({score,a1.residue,a2.residue});
            }
        }
        std::sort(pair_classes.begin(),pair_classes.end(),[](const PairClass&a,const PairClass&b){
            if (a.base_score!=b.base_score) return a.base_score<b.base_score;
            if (a.r1!=b.r1) return a.r1<b.r1;
            return a.r2<b.r2;
        });
        if ((int)pair_classes.size()>class_cap) pair_classes.resize(class_cap);
        std::cerr << "pair " << q1_fixed << ',' << q2 << " classes=" << pair_classes.size();
        if (!pair_classes.empty()) std::cerr << " best=" << pair_classes[0].base_score;
        std::cerr << "\n";

        for (const PairClass& cl:pair_classes) {
            ++classes;
            std::vector<int> rs=original;
            rs[q1_fixed]=cl.r1;
            rs[q2]=cl.r2;
            auto [x,M]=crt(mods,rs);
            u128 t0=x>=interval_lo ? 0 : ceil_div(interval_lo-x,M);
            for (u128 n=x+t0*M; n<=interval_hi; n+=M) {
                ++evaluated;
                Mask live=ALL;
                for (int i=0; i<NS; ++i) if (shifts[i]>n) clrbit(live,i);
                int r8=(int)(n%8);
                Mask m8{};
                for (int i=0; i<NS; ++i) {
                    if (!((live[i>>6]>>(i&63))&1ULL)) continue;
                    int y=(r8+8-(int)(shifts[i]%8))%8;
                    if (y==3 || y==6 || y==7) setbit(m8,i);
                }
                mask_andnot(live,m8);
                for (const CT& f:filters200) mask_andnot(live,f.cov[(u64)(n%f.mod)]);
                int score=popcount(live);
                if (score>48) continue;
                ++passed48;
                kept.push_back({score,0,cl.base_score,q1_fixed,cl.r1,q2,cl.r2,n,live});
                if ((int)kept.size()>keep200*8) {
                    std::nth_element(kept.begin(),kept.begin()+keep200,kept.end(),[](const Candidate&a,const Candidate&b){return worse_candidate(b,a);});
                    kept.resize(keep200);
                }
            }
        }
    }

    std::sort(kept.begin(),kept.end(),[](const Candidate&a,const Candidate&b){return worse_candidate(b,a);});
    if ((int)kept.size()>keep200) kept.resize(keep200);
    std::sort(kept.begin(),kept.end(),[](const Candidate&a,const Candidate&b){
        if (a.score200!=b.score200) return a.score200<b.score200;
        return a.n<b.n;
    });
    kept.erase(std::unique(kept.begin(),kept.end(),[](const Candidate&a,const Candidate&b){return a.n==b.n;}),kept.end());

    std::vector<int> refine_primes=primes_3mod4(200,1000000);
    #pragma omp parallel for schedule(dynamic,1)
    for (long long ci=0; ci<(long long)kept.size(); ++ci) {
        Candidate& c=kept[(std::size_t)ci];
        for (int w=0; w<WORDS; ++w) {
            u64 bits=c.live[w];
            while (bits) {
                int bit=__builtin_ctzll(bits);
                bits&=bits-1;
                int idx=64*w+bit;
                if (idx>=NS || shifts[idx]>c.n) { clrbit(c.live,idx); continue; }
                u128 r=c.n-shifts[idx];
                bool obstructed=false;
                for (int p:refine_primes) {
                    if (r%(u64)p) continue;
                    int e=0;
                    do { r/=(u64)p; ++e; } while (r%(u64)p==0);
                    if (e&1) { obstructed=true; break; }
                }
                if (obstructed) clrbit(c.live,idx);
            }
        }
        c.score1m=popcount(c.live);
    }

    std::sort(kept.begin(),kept.end(),[](const Candidate&a,const Candidate&b){
        if (a.score1m!=b.score1m) return a.score1m<b.score1m;
        if (a.score200!=b.score200) return a.score200<b.score200;
        return a.n<b.n;
    });
    if ((int)kept.size()>out1m) kept.resize(out1m);

    std::ofstream out(output);
    for (const Candidate& c:kept) {
        out << c.score1m << ' ' << c.score200 << ' ' << c.base_score << ' '
            << c.q1 << ' ' << c.r1 << ' ' << c.q2 << ' ' << c.r2 << ' '
            << str128(c.n);
        for (u64 word:c.live) out << ' ' << word;
        out << '\n';
    }

    std::cerr << "DONE q1=" << q1_fixed << " classes=" << classes
              << " evaluated=" << evaluated << " p200_pass=" << passed48
              << " retained=" << kept.size();
    if (!kept.empty()) std::cerr << " best1m=" << kept[0].score1m;
    std::cerr << "\n";
    return 0;
}
