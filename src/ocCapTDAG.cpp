#include "ocCapTDAG.h"
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <cassert>

using std::vector;
using std::string;

// -------- logging --------
static std::ofstream g_log;
static int g_vl = 2;
template<typename F> static inline void LOG(int lvl, F&& f){
    if (g_log.is_open() && g_vl>=lvl){ f(g_log); g_log.flush(); }
}
static inline void LOGOPEN(const string& path, bool append, int v){
    if (!g_log.is_open()){
        g_log.open(path, append?(std::ios::out|std::ios::app):(std::ios::out|std::ios::trunc));
        g_vl = v;
    }
}

// -------- small utils --------
static inline long long llabsll(long long x){ return x>=0?x:-x; }

// 128-bit key + hash
using u128 = unsigned __int128;
struct U128Hash {
    size_t operator()(const u128& x) const noexcept {
        uint64_t lo = (uint64_t)x;
        uint64_t hi = (uint64_t)(x >> 64);
        uint64_t h = lo ^ (hi * 0x9e3779b97f4a7c15ull);
        return (size_t)h;
    }
};

// encode helpers on base-B digits
static inline int digit_at(u128 key, int idx, const vector<u128>& powB, u128 B){
    // powB[idx] = B^(m-1-idx)
    u128 q = key / powB[idx];
    u128 r = q % B;
    return (int) (uint64_t) r;
}
static inline u128 enc_inc(u128 key, int idx, const vector<u128>& powB){
    return key + powB[idx];
}

// ====== ocCapTDAG impl ======

long long ocCapTDAG::weight_ij(int i,int j,int m,int h){
    long long w=0;
    if(i==0)   w -= 1;
    if(i==m-1) w += 1;
    if(j==0)   w -= 1;
    if(j==h-1) w += 1;
    return w;
}
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

// windows closure with Δ≤T
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
        // 子受父（正向）
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
        // 父受子（反向）
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

CapTDAGResult ocCapTDAG::solve(int m, int h){
    LOGOPEN(cfg_.logfile, cfg_.log_append, cfg_.vlevel);

    const int n = m*h;
    int T = cfg_.T>0? cfg_.T : n;

    // base and powB (128-bit safe)
    const u128 B = (u128)(h + 1);
    vector<u128> powB(m);
    powB[m-1] = (u128)1;
    for(int k=m-2;k>=0;--k){
        __uint128_t t = powB[k+1] * B;
        powB[k] = (u128)t;
    }

    LOG(1, [&](std::ofstream& o){
        o<<"[CapT-DAG/APT] m="<<m<<" h="<<h<<" n="<<n<<" T="<<T
         <<" v="<<cfg_.vlevel<<" K="<<cfg_.max_labels_per_key<<"\n";
    });

    // windows
    vector<vector<int>> LB, UB;
    compute_windows_spanT(m,h,T,n,LB,UB);
    LOG(2, [&](std::ofstream& o){ o<<"[SpanT] window sanity check passed\n"; });

    // feasibility lower bound from windows (necessary)
    auto compute_Tmin_required = [&](const vector<vector<int>>& LBx,
                                     const vector<vector<int>>& UBx)->int{
        int Tmin=1;
        for(int i=0;i<m;++i) for(int j=0;j+1<h;++j){
            Tmin = std::max(Tmin, LBx[i][j+1]-UBx[i][j]);
        }
        for(int i=0;i+1<m;++i) for(int j=0;j<h;++j){
            Tmin = std::max(Tmin, LBx[i+1][j]-UBx[i][j]);
        }
        return std::max(Tmin, 1);
    };
    int Tmin_req = compute_Tmin_required(LB,UB);
    LOG(1, [&](std::ofstream& o){ o<<"[Feas] Tmin_required="<<Tmin_req<<"\n"; });
    if(T < Tmin_req){
        LOG(1, [&](std::ofstream& o){
            o<<"[Feas] Given T="<<T<<" < Tmin_required => infeasible under OC\n";
        });
        throw std::runtime_error("Infeasible: T below OC-required minimum");
    }

    // Label
    struct Label {
        long long dist=0;
        u128 prev_key=0;
        int prev_row=-1;
        int prev_label_idx=-1;
        std::vector<uint16_t> rlast; // size m
        std::vector<uint16_t> clast; // size h
    };
    auto ensure_dims = [&](Label& L){
        if((int)L.rlast.size()!=m) L.rlast.assign(m, 0);
        if((int)L.clast.size()!=h) L.clast.assign(h, 0);
    };

    using LVec = vector<Label>;
    std::unordered_map<u128, LVec, U128Hash> cur, nxt;

    auto dominates = [&](const Label& B, const Label& A)->bool{
        if (B.rlast.size()!=A.rlast.size() || B.clast.size()!=A.clast.size()) return false;
        if (B.dist > A.dist) return false;
        bool strictly = (B.dist < A.dist);
        for (size_t i=0;i<A.rlast.size();++i){
            if (B.rlast[i] > A.rlast[i]) return false;   // small is better
            if (B.rlast[i] < A.rlast[i]) strictly = true;
        }
        for (size_t j=0;j<A.clast.size();++j){
            if (B.clast[j] > A.clast[j]) return false;
            if (B.clast[j] < A.clast[j]) strictly = true;
        }
        return strictly;
    };

    auto insert_label = [&](LVec& vec, Label&& L, size_t& pruned_dom, size_t& truncatedK){
        ensure_dims(L);
        // dominated by existing?
        for (const auto& e : vec) if (dominates(e, L)) { ++pruned_dom; return; }
        // rebuild out set removing those dominated by L
        LVec out; out.reserve(vec.size()+1);
        for (auto &e : vec){
            if (dominates(L, e)) { ++pruned_dom; continue; }
            ensure_dims(e);
            out.push_back(std::move(e));
        }
        out.push_back(std::move(L));

        int K = cfg_.max_labels_per_key;
        if (K>0 && (int)out.size()>K){
            std::nth_element(out.begin(), out.begin()+K, out.end(),
                            [](const Label& a, const Label& b){ return a.dist < b.dist; });
            out.resize(K);
            ++truncatedK;
        }
        vec.swap(out);
    };

    // init
    u128 key0 = (u128)0;
    {
        Label L0;
        L0.dist=0; L0.prev_key=0; L0.prev_row=-1; L0.prev_label_idx=-1;
        L0.rlast.assign(m, 0);
        L0.clast.assign(h, 0);
        cur.emplace(key0, LVec{std::move(L0)});
    }

    // necessary prune: row/col crowding (EDD necessary condition)
    auto row_crowding_prune = [&](const Label& lab, const vector<int>& a_prime, int Lnext)->bool{
        vector<int> d; d.reserve(m);
        for(int i=0;i<m;++i){
            if(a_prime[i] < h){
                int di = (int)lab.rlast[i] + T; // deadline for next row service
                d.push_back(di);
            }
        }
        if(d.empty()) return false;
        std::sort(d.begin(), d.end());
        for(int k=1;k<= (int)d.size(); ++k){
            if(d[k-1] < Lnext + k) return true; // impossible to serve k rows in next k steps
        }
        return false;
    };
    auto col_crowding_prune = [&](const Label& lab, const vector<int>& a_prime, int Lnext)->bool{
        vector<int> d; d.reserve(h);
        // placed count per column: #rows with a'[i] > j
        for(int j=0;j<h;++j){
            int placed = 0;
            for(int i=0;i<m;++i) if(a_prime[i] > j) ++placed;
            if(placed < m){
                int dj = (int)lab.clast[j] + T; // deadline for next col service
                d.push_back(dj);
            }
        }
        if(d.empty()) return false;
        std::sort(d.begin(), d.end());
        for(int k=1;k<= (int)d.size(); ++k){
            if(d[k-1] < Lnext + k) return true;
        }
        return false;
    };

    // keep layers for reconstruction
    vector< std::unordered_map<u128, LVec, U128Hash> > LV(n+1);

    for(int L=0; L<n; ++L){
        nxt.clear();
        nxt.reserve(std::max<size_t>(cur.size()*2, 16));

        size_t cand=0, blk_oc=0, blk_win=0, blk_tL=0, blk_tU=0, blk_dead=0, accepted=0;
        size_t pruned_dom=0, truncatedK=0;

        for(auto &kv : cur){
            u128 key = kv.first;
            const LVec& labels = kv.second;

            vector<int> a(m);
            for(int i=0;i<m;++i) a[i]=digit_at(key,i,powB,B);

            for(size_t li=0; li<labels.size(); ++li){
                const Label& lab = labels[li];

                for(int i=0;i<m;++i){
                    int j = a[i];
                    if(j>=h) continue;
                    ++cand;

                    // OC（非增）：a[i] + 1 <= a[i-1]
                    if(i>0 && j+1> a[i-1]){ ++blk_oc; continue; }

                    int Lnext = L+1;
                    // window
                    if(!(LB[i][j] <= Lnext && Lnext <= UB[i][j])){ ++blk_win; continue; }

                    // Δ≤T：左/上父
                    if(j>0){
                        int yL = (int)lab.rlast[i];
                        int dL = Lnext - yL;
                        if(dL<1 || dL>T){ ++blk_tL; continue; }
                    }
                    if(i>0){
                        int yU = (int)lab.clast[j];
                        int dU = Lnext - yU;
                        if(dU<1 || dU>T){ ++blk_tU; continue; }
                    }

                    // candidate
                    u128 key2 = enc_inc(key,i,powB);
                    long long step = (long long)Lnext * weight_ij(i,j,m,h) * cfg_.dV;
                    long long d1 = lab.dist + step;

                    Label Lnew;
                    Lnew.dist=d1; Lnew.prev_key=key; Lnew.prev_row=i; Lnew.prev_label_idx=(int)li;
                    Lnew.rlast = lab.rlast; Lnew.clast = lab.clast;
                    if((int)Lnew.rlast.size()!=m) Lnew.rlast.assign(m,0);
                    if((int)Lnew.clast.size()!=h) Lnew.clast.assign(h,0);
                    Lnew.rlast[i] = (uint16_t)Lnext;
                    Lnew.clast[j] = (uint16_t)Lnext;

                    // EDD 必要剪枝（安全）
                    vector<int> a_prime = a; a_prime[i] = j+1;
                    if (row_crowding_prune(Lnew, a_prime, Lnext) ||
                        col_crowding_prune(Lnew, a_prime, Lnext)) {
                        ++blk_dead; continue;
                    }

                    auto &vec = nxt[key2];
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
                o<<"[DP-stat] L="<<L
                 <<" cand="<<cand
                 <<" blk_oc="<<blk_oc<<" blk_win="<<blk_win
                 <<" blk_tL="<<blk_tL<<" blk_tU="<<blk_tU
                 <<" blk_dead="<<blk_dead
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

        LV[L] = std::move(cur);
        cur.swap(nxt);
    }
    LV[n] = cur;

    // terminal key
    u128 keyN=(u128)0; for(int i=0;i<m;++i) keyN += powB[i]*(u128)h;
    auto it = LV[n].find(keyN);
    if(it==LV[n].end() || it->second.empty()){
        LOG(1, [&](std::ofstream& o){ o<<"[END] terminal state missing\n"; });
        throw std::runtime_error("No feasible terminal state");
    }
    int best_idx=0; for(int t=1;t<(int)it->second.size();++t) if(it->second[t].dist < it->second[best_idx].dist) best_idx=t;

    // reconstruct
    vector<vector<int>> y(m, vector<int>(h,0));
    u128 k = keyN; int li = best_idx;
    for(int L=n; L>=1; --L){
        const auto &mp = LV[L];
        auto hit = mp.find(k);
        if(hit==mp.end()) throw std::runtime_error("Reconstruct failed (key missing)");
        const vector<Label>& V = hit->second;
        if(li<0 || li>=(int)V.size()) li = 0;
        const Label& curLab = V[li];

        u128 pk = curLab.prev_key;
        int ri = curLab.prev_row;
        int aj = digit_at(pk, ri, powB, B);
        if(ri<0 || ri>=m || aj<0 || aj>=h) throw std::runtime_error("Reconstruct index OOR");
        y[ri][aj] = L;

        long long step_cost = (long long)L * weight_ij(ri, aj, m, h) * cfg_.dV;

        int next_li = curLab.prev_label_idx;
        if(L-1>=0){
            const auto &mp2 = LV[L-1];
            auto it2 = mp2.find(pk);
            if(it2==mp2.end() || it2->second.empty()){
                throw std::runtime_error("Reconstruct failed (prev state missing)");
            }
            const vector<Label>& PV = it2->second;
            if(next_li<0 || next_li>=(int)PV.size()){
                int found=-1;
                for(int t=0;t<(int)PV.size();++t){
                    if(PV[t].dist + step_cost == curLab.dist){ found = t; break; }
                }
                if(found<0){
                    long long best = PV[0].dist;
                    found = 0;
                    for(int t=1;t<(int)PV.size();++t){
                        if(PV[t].dist < best){ best=PV[t].dist; found=t; }
                    }
                }
                next_li = found;
            }
        }
        k = curLab.prev_key;
        li = next_li;
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
                for(int j=0;j<h;++j){
                    o<< (j?" ":"") << (y[i][j]<10?"  ":" ") << y[i][j];
                }
                o<<"\n";
            }
        });
    }
    if(T>0 && dmax > T){
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
