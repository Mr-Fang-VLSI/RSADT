#include "ocCapTDAG.h"
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <limits>
#include <cassert>

using std::vector;
static std::ofstream g_log;
static int g_vl=2;
static inline void LOGOPEN(const std::string& path, bool append, int v){ if(!g_log.is_open()){ g_log.open(path, append?(std::ios::out|std::ios::app):(std::ios::out|std::ios::trunc)); g_vl=v; } }
template<typename F> static inline void LOG(int lvl, F&& f){ if(g_log.is_open() && g_vl>=lvl){ f(g_log); g_log.flush(); } }

long long ocCapTDAG::weight_ij(int i,int j,int m,int h){
    long long w=0;
    if(i==0)   w -= 1;
    if(i==m-1) w += 1;
    if(j==0)   w -= 1;
    if(j==h-1) w += 1;
    return w;
}
static inline long long llabsll(long long x){ return x>=0?x:-x; }

long long ocCapTDAG::hpwl_neighbors(const vector<vector<int>>& y, long long dV){
    const int m=(int)y.size(), h=(int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) S += llabsll((long long)y[i][j+1]-y[i][j])*dV;
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j)   S += llabsll((long long)y[i+1][j]-y[i][j])*dV;
    return S;
}
bool ocCapTDAG::check_OC(const vector<vector<int>>& y){
    const int m=(int)y.size(), h=(int)y[0].size();
    for(int i1=0;i1<m;++i1) for(int j1=0;j1<h;++j1){
        int v1=y[i1][j1]; if(v1<=0) return false;
        for(int i2=i1;i2<m;++i2) for(int j2=j1;j2<h;++j2) if(y[i1][j1]>y[i2][j2]) return false;
    }
    return true;
}
int ocCapTDAG::max_delta_adj(const vector<vector<int>>& y){
    const int m=(int)y.size(), h=(int)y[0].size();
    int dmax=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) dmax=std::max(dmax, y[i][j+1]-y[i][j]);
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j) dmax=std::max(dmax, y[i+1][j]-y[i][j]);
    return dmax;
}

// Δ≤T necessary + OC closure on [LB,UB]
void ocCapTDAG::compute_windows_spanT(int m,int h,int T,int n,
                                      vector<vector<int>>& LB,
                                      vector<vector<int>>& UB) const
{
    LB.assign(m, vector<int>(h,0));
    UB.assign(m, vector<int>(h,0));
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        LB[i][j] = (i+1)*(j+1);
        UB[i][j] = n - (m-1-i)*(h-1-j);
    }
    if(T<=0) T=n;

    bool changed=true; int iter=0, iter_max=m*h*6;
    while(changed && ++iter<=iter_max){
        changed=false;
        // 子受父
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            int oldL=LB[i][j], oldU=UB[i][j];
            if(i>0)   LB[i][j]=std::max(LB[i][j], LB[i-1][j]+1);
            if(j>0)   LB[i][j]=std::max(LB[i][j], LB[i][j-1]+1);
            if(i>0)   UB[i][j]=std::min(UB[i][j], UB[i-1][j]+T);
            if(j>0)   UB[i][j]=std::min(UB[i][j], UB[i][j-1]+T);
            if(i<m-1) UB[i][j]=std::min(UB[i][j], UB[i+1][j]-1);
            if(j<h-1) UB[i][j]=std::min(UB[i][j], UB[i][j+1]-1);
            if(LB[i][j]!=oldL || UB[i][j]!=oldU) changed=true;
        }
        // 父受子
        for(int i=m-1;i>=0;--i) for(int j=h-1;j>=0;--j){
            int oldL=LB[i][j], oldU=UB[i][j];
            if(i<m-1) LB[i][j]=std::max(LB[i][j], LB[i+1][j]-T);
            if(j<h-1) LB[i][j]=std::max(LB[i][j], LB[i][j+1]-T);
            if(i<m-1) UB[i][j]=std::min(UB[i][j], UB[i+1][j]-1);
            if(j<h-1) UB[i][j]=std::min(UB[i][j], UB[i][j+1]-1);
            if(LB[i][j]!=oldL || UB[i][j]!=oldU) changed=true;
        }
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            if(LB[i][j]>UB[i][j]){
                LOG(1, [&](std::ofstream& o){ o<<"[SpanT] infeasible at ("<<i<<","<<j<<") LB="<<LB[i][j]<<" UB="<<UB[i][j]<<"\n"; });
                throw std::runtime_error("Δ≤T infeasible after closure");
            }
        }
    }
    if(cfg_.vlevel>=3){
        int LBmin=1e9,LBmax=-1e9,UBmin=1e9,UBmax=-1e9;
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            LBmin=std::min(LBmin,LB[i][j]); LBmax=std::max(LBmax,LB[i][j]);
            UBmin=std::min(UBmin,UB[i][j]); UBmax=std::max(UBmax,UB[i][j]);
        }
        LOG(3, [&](std::ofstream& o){
            o<<"[SpanT] LB range=["<<LBmin<<","<<LBmax<<"], UB range=["<<UBmin<<","<<UBmax<<"]\n";
        });
    }
}

// -------- solve: 单向前推 + cap=0 边 ----------
CapTDAGResult ocCapTDAG::solve(int m, int h){
    LOGOPEN(cfg_.logfile, cfg_.log_append, cfg_.vlevel);

    const int n = m*h;
    int T = cfg_.T>0? cfg_.T : n;
    const uint64_t B = (uint64_t)h + 1u;

    // precompute powB
    vector<uint64_t> powB(m);
    powB[m-1]=1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t t=(__uint128_t)powB[k+1]*(__uint128_t)B;
        if(t > (__uint128_t)std::numeric_limits<uint64_t>::max()) throw std::runtime_error("state key overflow");
        powB[k]=(uint64_t)t;
    }

    LOG(1, [&](std::ofstream& o){
        o<<"[CapT-DAG] m="<<m<<" h="<<h<<" n="<<n<<" T="<<T<<" v="<<cfg_.vlevel<<"\n";
    });

    // windows
    vector<vector<int>> LB, UB;
    compute_windows_spanT(m,h,T,n,LB,UB);
    LOG(2, [&](std::ofstream& o){ o<<"[SpanT] window sanity check passed\n"; });

    // DP maps
    struct Node { long long dist; uint64_t prev; int prev_row; };
    std::unordered_map<uint64_t, Node> cur, nxt;
    cur.reserve(1024); nxt.reserve(4096);
    vector< std::unordered_map<uint64_t, std::pair<uint64_t,int>> > parents(n+1);
    for(auto& mp : parents) mp.reserve(1024);

    uint64_t key0=0ull;
    cur.emplace(key0, Node{0,0ull,-1});

    auto valid_trans = [&](int i,int j,int L)->bool{
        if(!(LB[i][j] <= L+1 && L+1 <= UB[i][j])) return false;
        if(j>0){ if(L+1 > LB[i][j-1] + T) return false; } // Δ upper bound via parent-LB
        if(i>0){ if(L+1 > LB[i-1][j] + T) return false; }
        return true;
    };

    for(int L=0; L<n; ++L){
        nxt.clear();
        nxt.reserve(std::max<size_t>(cur.size()*2, 16));
        for(auto &kv : cur){
            uint64_t key = kv.first;
            long long d0 = kv.second.dist;
            // decode skyline
            int a[64]; for(int i=0;i<m;++i) a[i]=digit_at(key,i,powB,B);
            for(int i=0;i<m;++i){
                int j = a[i];
                if(j>=h) continue;
                if(i>0 && j+1> a[i-1]) continue; // OC
                if(!valid_trans(i,j,L)) continue; // cap=0 边: 直接丢弃

                uint64_t key2 = enc_inc(key,i,powB);
                long long cost = (long long)(L+1) * weight_ij(i,j,m,h) * cfg_.dV;
                long long d1 = d0 + cost;

                auto it = nxt.find(key2);
                if(it==nxt.end()){
                    nxt.emplace(key2, Node{d1, key, i});
                    parents[L+1].emplace(key2, std::make_pair(key, i));
                }else if(d1 < it->second.dist){
                    it->second.dist = d1;
                    it->second.prev = key;
                    it->second.prev_row = i;
                    parents[L+1][key2] = std::make_pair(key, i);
                }
            }
        }
        if(cfg_.progress && cfg_.vlevel>=2 && (L%16==0 || L+1==n)){
            LOG(2, [&](std::ofstream& o){
                o<<"[DP] level "<<L<<" states="<<cur.size()<<" -> next="<<nxt.size()<<"\n";
            });
        }
        if(nxt.empty()){
            LOG(1, [&](std::ofstream& o){ o<<"[DP] next empty at level "<<L<<"\n"; });
            throw std::runtime_error("No feasible next frontier");
        }
        cur.swap(nxt);
    }

    // terminal state
    uint64_t keyN=0ull; for(int i=0;i<m;++i) keyN += powB[i]*(uint64_t)h;
    auto it = cur.find(keyN);
    if(it==cur.end()){
        LOG(1, [&](std::ofstream& o){ o<<"[END] terminal state missing\n"; });
        throw std::runtime_error("No feasible terminal state");
    }

    // reconstruct path
    vector<vector<int>> y(m, vector<int>(h,0));
    uint64_t k = keyN;
    for(int L=n; L>=1; --L){
        auto pit = parents[L].find(k);
        if(pit==parents[L].end()) throw std::runtime_error("Reconstruct failed");
        uint64_t pk = pit->second.first;
        int ri = pit->second.second;
        // j = a_i before inc
        int aj = digit_at(pk, ri, powB, B);
        y[ri][aj] = L;
        k = pk;
    }

    long long total = hpwl_neighbors(y, cfg_.dV);
    int dmax = max_delta_adj(y);
    bool oc_ok = check_OC(y);

    LOG(1, [&](std::ofstream& o){
        o<<"[Verify] HPWL="<<total<<" maxΔ="<<dmax<<" OC="<<(oc_ok?"OK":"FAIL")<<"\n";
    });
    if(cfg_.vlevel>=3){
        LOG(3, [&](std::ofstream& o){
            o<<"[y_order top->bottom]\n";
            for(int i=0;i<m;++i){
                o<<"  ";
                for(int j=0;j<h;++j) o<< (j?" ":"") << (y[i][j]<10?"  ":" ") << y[i][j];
                o<<"\n";
            }
        });
    }

    if(dmax > T){
        LOG(1, [&](std::ofstream& o){ o<<"[POST] maxΔ="<<dmax<<" > T="<<T<<"\n"; });
        throw std::runtime_error("Post-check failed: max Δ > T (should not happen with cap-edges)");
    }

    CapTDAGResult R;
    R.m=m; R.h=h; R.n=n; R.T=T; R.dV=cfg_.dV;
    R.y_order=std::move(y);
    R.hpwl=total;
    R.max_delta=dmax;
    R.oc_ok=oc_ok;
    return R;
}
