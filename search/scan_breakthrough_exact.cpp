#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_set>
#include <vector>

using u64 = std::uint64_t;
using u128 = __uint128_t;

static inline u64 mul_mod(u64 a,u64 b,u64 m){return static_cast<u64>((static_cast<u128>(a)*b)%m);}
static inline u64 pow_mod(u64 a,u64 e,u64 m){u64 r=1;while(e){if(e&1U)r=mul_mod(r,a,m);a=mul_mod(a,a,m);e>>=1U;}return r;}
static bool is_prime_u64(u64 n){
 if(n<2)return false;
 for(u64 p:{2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL,29ULL,31ULL,37ULL})if(n%p==0)return n==p;
 u64 d=n-1;unsigned s=0;while((d&1U)==0){d>>=1U;++s;}
 for(u64 a:{2ULL,325ULL,9375ULL,28178ULL,450775ULL,9780504ULL,1795265022ULL}){
  if(a%n==0)continue;u64 x=pow_mod(a%n,d,n);if(x==1||x==n-1)continue;bool witness=true;
  for(unsigned r=1;r<s;++r){x=mul_mod(x,x,n);if(x==n-1){witness=false;break;}}
  if(witness)return false;
 }
 return true;
}
static inline u64 splitmix64(u64&state){u64 z=(state+=0x9e3779b97f4a7c15ULL);z=(z^(z>>30U))*0xbf58476d1ce4e5b9ULL;z=(z^(z>>27U))*0x94d049bb133111ebULL;return z^(z>>31U);}
static u64 brent_factor(u64 n){
 if((n&1U)==0)return 2;if(n%3==0)return 3;
 for(u64 attempt=0;;++attempt){
  u64 state=n^(attempt*0xd6e8feb86659fd93ULL);u64 y=splitmix64(state)%(n-1)+1,c=splitmix64(state)%(n-1)+1,m=128,g=1,r=1,q=1,x=0,ys=0,work=0;constexpr u64 LIMIT=1ULL<<21U;
  auto f=[&](u64 v){return(mul_mod(v,v,n)+c)%n;};
  while(g==1&&work<LIMIT){x=y;for(u64 i=0;i<r&&work<LIMIT;++i,++work)y=f(y);for(u64 k=0;k<r&&g==1&&work<LIMIT;k+=m){ys=y;q=1;u64 lim=std::min(m,r-k);for(u64 i=0;i<lim&&work<LIMIT;++i,++work){y=f(y);u64 diff=x>y?x-y:y-x;q=mul_mod(q,diff==0?n:diff,n);}g=std::gcd(q,n);}r<<=1U;}
  if(work>=LIMIT||g==1)continue;
  if(g==n){u64 fallback=0;do{ys=f(ys);u64 diff=x>ys?x-ys:ys-x;g=std::gcd(diff,n);++fallback;}while(g==1&&fallback<LIMIT);if(g==n||g==1||fallback>=LIMIT)continue;}
  return g;
 }
}
static void factor_u64(u64 n,std::vector<u64>&out){if(n==1)return;if(is_prime_u64(n)){out.push_back(n);return;}u64 d=brent_factor(n);factor_u64(d,out);factor_u64(n/d,out);}
static std::vector<int> make_small_primes(int limit){std::vector<bool>ok(limit+1,true);ok[0]=ok[1]=false;std::vector<int>ps;for(int p=2;p<=limit;++p)if(ok[p]){ps.push_back(p);if(1LL*p*p<=limit)for(int q=p*p;q<=limit;q+=p)ok[q]=false;}return ps;}
static const std::vector<int> SMALL_PRIMES=make_small_primes(10000);
static bool is_sum_two_squares(u64 n){
 if(n==0)return true;
 for(int p:SMALL_PRIMES){if(static_cast<u64>(p)*p>n)break;if(n%static_cast<u64>(p))continue;unsigned e=0;do{n/=static_cast<u64>(p);++e;}while(n%static_cast<u64>(p)==0);if((p&3)==3&&(e&1U))return false;}
 if(n==1)return true;if(is_prime_u64(n))return(n&3U)!=3U;
 std::vector<u64>fs;fs.reserve(16);factor_u64(n,fs);std::sort(fs.begin(),fs.end());
 for(std::size_t i=0;i<fs.size();){std::size_t j=i+1;while(j<fs.size()&&fs[j]==fs[i])++j;if((fs[i]&3U)==3U&&((j-i)&1U))return false;i=j;}return true;
}
static inline u64 isqrt_u64(u64 n){u64 x=static_cast<u64>(std::sqrt(static_cast<long double>(n)));while(static_cast<u128>(x+1)*(x+1)<=n)++x;while(static_cast<u128>(x)*x>n)--x;return x;}
static inline u64 ceil_sqrt_u64(u64 n){u64 r=isqrt_u64(n);return static_cast<u128>(r)*r==n?r:r+1;}
static std::vector<u64> make_B(u64 limit){std::vector<u64>B{1};for(int k=1;;++k){u128 next=static_cast<u128>(B.back())*(2*k)*(2*k+1)/(static_cast<u64>(k)*(k+1));if(next>limit)break;B.push_back(static_cast<u64>(next));}return B;}
static std::vector<u64> make_shifts_for_K(int K){std::vector<u64>B(K+1);B[0]=1;for(int k=1;k<=K;++k)B[k]=static_cast<u64>(static_cast<u128>(B[k-1])*(2*k)*(2*k+1)/(static_cast<u64>(k)*(k+1)));std::vector<u64>S;for(int c=0;c<=K;++c)for(int d=c;d<=K;++d)S.push_back(B[c]+B[d]);std::sort(S.begin(),S.end());S.erase(std::unique(S.begin(),S.end()),S.end());return S;}
static std::vector<u64> make_all_shifts(u64 limit){auto B=make_B(limit);std::vector<u64>S;for(std::size_t c=0;c<B.size();++c)for(std::size_t d=c;d<B.size();++d){u128 s=static_cast<u128>(B[c])+B[d];if(s<=limit)S.push_back(static_cast<u64>(s));}std::sort(S.begin(),S.end());S.erase(std::unique(S.begin(),S.end()),S.end());return S;}
static inline void set_bit(std::vector<u64>&b,u64 i){b[i>>6]|=1ULL<<(i&63);}static inline bool get_bit(const std::vector<u64>&b,u64 i){return(b[i>>6]>>(i&63))&1ULL;}
static inline u64 extract_low(const std::vector<u64>&src,u64 pos,unsigned take){std::size_t w=static_cast<std::size_t>(pos>>6);unsigned off=pos&63;u64 x=src[w]>>off;if(off&&w+1<src.size())x|=src[w+1]<<(64-off);if(take<64)x&=(1ULL<<take)-1;return x;}
static void or_bit_range(std::vector<u64>&dst,u64 dpos,const std::vector<u64>&src,u64 spos,u64 len){while(len){unsigned doff=dpos&63;unsigned take=static_cast<unsigned>(std::min<u64>(len,64-doff));u64 x=extract_low(src,spos,take);dst[dpos>>6]|=x<<doff;dpos+=take;spos+=take;len-=take;}}
static std::vector<u64> exact_s2_interval_lattice(u64 lo,u64 hi,u64&pair_visits){
 const u64 len=hi-lo+1;std::vector<u64>good(static_cast<std::size_t>((len+63)>>6),0);u64 bhi=isqrt_u64(hi),blo=ceil_sqrt_u64(lo),amax=isqrt_u64(hi/2);
 for(u64 a=0;a<=amax;++a){const u128 aa=static_cast<u128>(a)*a;while(bhi>=a&&aa+static_cast<u128>(bhi)*bhi>hi)--bhi;if(blo<a)blo=a;while(blo>a&&aa+static_cast<u128>(blo-1)*(blo-1)>=lo)--blo;while(blo<=bhi&&aa+static_cast<u128>(blo)*blo<lo)++blo;if(blo>bhi)continue;for(u64 b=blo;b<=bhi;++b){u64 value=static_cast<u64>(aa+static_cast<u128>(b)*b);set_bit(good,value-lo);++pair_visits;}}
 return good;
}
int main(int argc,char**argv){
 if(argc<4){std::cerr<<"usage: scan_breakthrough START END OUTPUT [K0=12]\n";return 2;}
 const u64 nlo=std::strtoull(argv[1],nullptr,10),nhi=std::strtoull(argv[2],nullptr,10);const std::string output_path=argv[3];const int K0=argc>4?std::atoi(argv[4]):12;
 if(nlo>nhi||K0<0){std::cerr<<"invalid range or K0\n";return 2;}
 const u64 count=nhi-nlo+1;std::vector<u64>B(K0+1);B[0]=1;for(int k=1;k<=K0;++k)B[k]=static_cast<u64>(static_cast<u128>(B[k-1])*(2*k)*(2*k+1)/(static_cast<u64>(k)*(k+1)));
 const u64 bmin=B.front(),bmax=B.back();if(nlo<2*bmax){std::cerr<<"START is too small for this K0\n";return 2;}
 const u64 tlo=nlo-2*bmax,thi=nhi-2*bmin,ulo=nlo-bmax,uhi=nhi-bmin,tlen=thi-tlo+1,ulen=uhi-ulo+1;
 auto t0=std::chrono::steady_clock::now();u64 lattice_pairs=0;auto T=exact_s2_interval_lattice(tlo,thi,lattice_pairs);auto t1=std::chrono::steady_clock::now();
 std::vector<u64>U(static_cast<std::size_t>((ulen+63)>>6),0);for(u64 b:B)or_bit_range(U,0,T,(ulo-b)-tlo,ulen);
 std::vector<u64>covered(static_cast<std::size_t>((count+63)>>6),0);for(u64 b:B)or_bit_range(covered,0,U,(nlo-b)-ulo,count);auto t2=std::chrono::steady_clock::now();
 std::vector<u64>candidates;for(u64 i=0;i<count;++i)if(!get_bit(covered,i))candidates.push_back(nlo+i);
 auto small_shifts=make_shifts_for_K(K0),all_shifts=make_all_shifts(nhi);std::unordered_set<u64>already_checked(small_shifts.begin(),small_shifts.end());
 std::vector<u64>counterexamples;u64 exact_calls=0;for(u64 n:candidates){bool represented=false;for(u64 shift:all_shifts){if(shift>n)break;if(already_checked.count(shift))continue;++exact_calls;if(is_sum_two_squares(n-shift)){represented=true;break;}}if(!represented)counterexamples.push_back(n);}auto t3=std::chrono::steady_clock::now();
 double lattice_seconds=std::chrono::duration<double>(t1-t0).count(),dilation_seconds=std::chrono::duration<double>(t2-t1).count(),exact_seconds=std::chrono::duration<double>(t3-t2).count(),total=lattice_seconds+dilation_seconds+exact_seconds;
 std::ofstream output(output_path);if(!output){std::cerr<<"could not open output\n";return 2;}
 output<<"{\n  \"schema\": \"a303639-breakthrough-exact-v1\",\n  \"start\": "<<nlo<<",\n  \"end\": "<<nhi<<",\n  \"processed\": "<<count<<",\n  \"prefilter_K\": "<<K0<<",\n  \"binomial_shifts_per_stage\": "<<B.size()<<",\n  \"equivalent_pair_shifts\": "<<small_shifts.size()<<",\n  \"s2_interval_length\": "<<tlen<<",\n  \"lattice_pair_visits\": "<<lattice_pairs<<",\n  \"prefilter_candidates\": "<<candidates.size()<<",\n  \"exact_calls_after_prefilter\": "<<exact_calls<<",\n  \"lattice_seconds\": "<<std::setprecision(17)<<lattice_seconds<<",\n  \"factorized_dilation_seconds\": "<<dilation_seconds<<",\n  \"exact_seconds\": "<<exact_seconds<<",\n  \"rate_numbers_per_second\": "<<static_cast<double>(count)/total<<",\n  \"counterexamples\": [";
 for(std::size_t i=0;i<counterexamples.size();++i){if(i)output<<", ";output<<counterexamples[i];}output<<"]\n}\n";
 std::cerr<<"processed="<<count<<" candidates="<<candidates.size()<<" exact_calls="<<exact_calls<<" counterexamples="<<counterexamples.size()<<" rate="<<static_cast<double>(count)/total<<"\n";return 0;
}
