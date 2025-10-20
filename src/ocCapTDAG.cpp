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
        for(int i2=i1;i2<m;++i2) for(int j2=j1;j2<h;++j2)
            if(y[i1][j1]>y[i2][j2]) return false;
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

// -------- APT 多标签前向 DP ----------
CapTDAGResult ocCapTDAG::solve(int m, int h){
    LOGOPEN(cfg_.logfile, cfg_.log_append, cfg_.vlevel);

    const int n = m*h;
    int T = cfg_.T>0? cfg_.T : n;
    const uint64_t B = (uint64_t)h + 1u;

    // powB
    vector<uint64_t> powB(m);
    powB[m-1]=1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t t=(__uint128_t)powB[k+1]*(__uint128_t)B;
        if(t > (__uint128_t)std::numeric_limits<uint64_t>::max()) throw std::runtime_error("state key overflow");
        powB[k]=(uint64_t)t;
    }

    LOG(1, [&](std::ofstream& o){
        o<<"[CapT-DAG/APT] m="<<m<<" h="<<h<<" n="<<n<<" T="<<T<<" v="<<cfg_.vlevel
         <<" K="<<cfg_.max_labels_per_key<<"\n";
    });

    // windows
    vector<vector<int>> LB, UB;
    compute_windows_spanT(m,h,T,n,LB,UB);
    LOG(2, [&](std::ofstream& o){ o<<"[SpanT] window sanity check passed\n"; });

    // Label struct
    struct Label {
        long long dist=0;
        uint64_t prev_key=0;
        int prev_row=-1;
        int prev_label_idx=-1;
        vector<short> rlast; // size m
        vector<short> clast; // size h
    };

    // map: key -> vector<labels>
    using LVec = vector<Label>;
    std::unordered_map<uint64_t, LVec> cur, nxt;

    auto dominated_by = [&](const Label& A, const Label& B)->bool{
        // B dominates A if: B.dist <= A.dist and B.rlast >= A.rlast (all) and B.clast >= A.clast (all)
        if(B.dist > A.dist) return false;
        for(size_t i=0;i<A.rlast.size();++i) if(B.rlast[i] < A.rlast[i]) return false;
        for(size_t j=0;j<A.clast.size();++j) if(B.clast[j] < A.clast[j]) return false;
        // strictly better in at least one dimension
        return (B.dist < A.dist);
    };

    auto insert_label = [&](LVec& vec, Label&& L, size_t& pruned_dom, size_t& truncatedK){
        // if dominated by any -> drop
        for(const auto& e: vec){
            if(dominated_by(L, e)){ ++pruned_dom; return; }
        }
        // remove labels dominated by L
        size_t w=0;
        for(size_t t=0;t<vec.size();++t){
            if(dominated_by(vec[t], L)) { ++pruned_dom; continue; }
            vec[w++]=std::move(vec[t]);
        }
        vec.resize(w);
        // push
        vec.push_back(std::move(L));
        // cap by K (drop worst dist)
        int K = cfg_.max_labels_per_key;
        if(K>0 && (int)vec.size()>K){
            // find worst dist
            size_t worst=0;
            for(size_t t=1;t<vec.size();++t){
                if(vec[t].dist > vec[worst].dist) worst=t;
            }
            if(worst < vec.size()-1) std::swap(vec[worst], vec.back());
            vec.pop_back(); ++truncatedK;
        }
    };

    // init
    uint64_t key0=0ull;
    {
        Label L0;
        L0.dist=0; L0.prev_key=0; L0.prev_row=-1; L0.prev_label_idx=-1;
        L0.rlast.assign(m, 0);
        L0.clast.assign(h, 0);
        cur.emplace(key0, LVec{std::move(L0)});
    }

    // keep all levels for reconstruction
    vector< std::unordered_map<uint64_t, LVec> > LV(n+1);

    for(int L=0; L<n; ++L){
        nxt.clear();
        nxt.reserve(std::max<size_t>(cur.size()*2, 16));

        size_t cand=0, blk_oc=0, blk_win=0, blk_tL=0, blk_tU=0, accepted=0;
        size_t pruned_dom=0, truncatedK=0;

        for(auto &kv : cur){
            uint64_t key = kv.first;
            const LVec& labels = kv.second;

            int a[64]; for(int i=0;i<m;++i) a[i]=digit_at(key,i,powB,B);

            for(size_t li=0; li<labels.size(); ++li){
                const Label& lab = labels[li];

                for(int i=0;i<m;++i){
                    int j = a[i];
                    if(j>=h) continue;
                    ++cand;

                    // OC
                    if(i>0 && j+1> a[i-1]){ ++blk_oc; continue; }

                    int Lnext = L+1;
                    // window
                    if(!(LB[i][j] <= Lnext && Lnext <= UB[i][j])){ ++blk_win; continue; }

                    // T-left
                    if(j>0){
                        int yL = (int)lab.rlast[i];
                        int dL = Lnext - yL;
                        if(dL<1 || dL>T){ ++blk_tL; continue; }
                    }
                    // T-up
                    if(i>0){
                        int yU = (int)lab.clast[j];
                        int dU = Lnext - yU;
                        if(dU<1 || dU>T){ ++blk_tU; continue; }
                    }

                    uint64_t key2 = enc_inc(key,i,powB);
                    long long cost = (long long)(Lnext) * weight_ij(i,j,m,h) * cfg_.dV;
                    long long d1 = lab.dist + cost;

                    Label Lnew;
                    Lnew.dist=d1; Lnew.prev_key=key; Lnew.prev_row=i; Lnew.prev_label_idx=(int)li;
                    Lnew.rlast = lab.rlast; Lnew.clast = lab.clast;
                    Lnew.rlast[i] = (short)Lnext;
                    Lnew.clast[j] = (short)Lnext;

                    auto &vec = nxt[key2]; // creates empty if not exist
                    insert_label(vec, std::move(Lnew), pruned_dom, truncatedK);
                    ++accepted;
                }
            }
        }

        if(cfg_.progress && (cfg_.vlevel>=2) && (L%16==0 || L+1==n)){
            LOG(2, [&](std::ofstream& o){
                o<<"[DP] level "<<L<<" states="<<cur.size()<<" -> next="<<nxt.size()<<"\n";
            });
        }
        if(cfg_.vlevel>=3){
            size_t next_labels=0;
            for(auto &p: nxt) next_labels += p.second.size();
            double avg = nxt.empty()?0.0: (double)next_labels / (double)nxt.size();
            LOG(3, [&](std::ofstream& o){
                o<<"[DP-stat] L="<<L<<" cand="<<cand
                 <<" blk_oc="<<blk_oc<<" blk_win="<<blk_win
                 <<" blk_tL="<<blk_tL<<" blk_tU="<<blk_tU
                 <<" accepted="<<accepted
                 <<" pruned_dom="<<pruned_dom
                 <<" truncK="<<truncatedK
                 <<" next_keys="<<nxt.size()
                 <<" avg_labels/key="<<avg
                 <<"\n";
            });
        }

        if(nxt.empty()){
            LOG(1, [&](std::ofstream& o){ o<<"[DP] next empty at level "<<L<<"\n"; });
            throw std::runtime_error("No feasible next frontier");
        }

        // 保存本层，供回溯
        LV[L] = std::move(cur);
        cur.swap(nxt);
    }
    LV[n] = cur;

    // terminal state
    uint64_t keyN=0ull; for(int i=0;i<m;++i) keyN += powB[i]*(uint64_t)h;
    auto it = LV[n].find(keyN);
    if(it==LV[n].end() || it->second.empty()){
        LOG(1, [&](std::ofstream& o){ o<<"[END] terminal state missing\n"; });
        throw std::runtime_error("No feasible terminal state");
    }
    // 选 dist 最小的 label
    int best_idx=0; for(int t=1;t<(int)it->second.size();++t) if(it->second[t].dist < it->second[best_idx].dist) best_idx=t;

    // reconstruct path
    vector<vector<int>> y(m, vector<int>(h,0));
    uint64_t k = keyN; int li = best_idx;
    for(int L=n; L>=1; --L){
        const auto &mp = LV[L];
        auto hit = mp.find(k);
        if(hit==mp.end()) throw std::runtime_error("Reconstruct failed (key missing)");
        const vector<Label>& V = hit->second;
        if(li<0 || li>=(int)V.size()) throw std::runtime_error("Reconstruct failed (label idx)");
        const Label& curLab = V[li];

        uint64_t pk = curLab.prev_key;
        int ri = curLab.prev_row;
        int aj = 0;
        // decode aj from pk
        {
            const int a_i = digit_at(pk, ri, powB, B);
            aj = a_i; // before increment
        }
        y[ri][aj] = L;

        // step back
        k = pk;
        li = curLab.prev_label_idx;
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
        throw std::runtime_error("Post-check failed: max Δ > T");
    }

    CapTDAGResult R;
    R.m=m; R.h=h; R.n=n; R.T=T; R.dV=cfg_.dV;
    R.y_order=std::move(y);
    R.hpwl=total;
    R.max_delta=dmax;
    R.oc_ok=oc_ok;
    return R;
}
