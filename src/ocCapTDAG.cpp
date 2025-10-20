#include "ocCapTDAG.h"
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <limits>
#include <cstdint>
#include <climits>

using std::vector;

// ------------- Logger -------------
static std::ofstream g_log;
static int g_vl = 2;
template<typename F>
static inline void LOG(int lvl, F&& f){
    if (g_log.is_open() && g_vl >= lvl) { f(g_log); g_log.flush(); }
}
static inline void LOGOPEN(const std::string& path, bool append, int v){
    if (!g_log.is_open()){
        g_log.open(path, append ? (std::ios::out|std::ios::app) : (std::ios::out|std::ios::trunc));
        g_vl = v;
    }
}

static inline long long llabsll(long long x){ return x>=0?x:-x; }
static inline uint16_t sat_inc(uint16_t v, uint16_t cap){ return (v<cap)? (uint16_t)(v+1) : cap; }

// ------------- Helpers -------------
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
    for(int i1=0;i1<m;++i1) for(int j1=0;j1<h;++j1)
        for(int i2=i1;i2<m;++i2) for(int j2=j1;j2<h;++j2)
            if(y[i1][j1] > y[i2][j2]) return false;
    return true;
}
int ocCapTDAG::max_delta_adj(const vector<vector<int>>& y){
    const int m=(int)y.size(), h=(int)y[0].size();
    int dmax=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) dmax=std::max(dmax, y[i][j+1]-y[i][j]);
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j) dmax=std::max(dmax, y[i+1][j]-y[i][j]);
    return dmax;
}

// ------------- Δ≤T 窗口闭包（必要 + OC 单调）-------------
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
        // 子受父（必要 + OC）
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            int oL=LB[i][j], oU=UB[i][j];
            if(i>0)   LB[i][j] = std::max(LB[i][j], LB[i-1][j] + 1);
            if(j>0)   LB[i][j] = std::max(LB[i][j], LB[i][j-1] + 1);
            if(i>0)   UB[i][j] = std::min(UB[i][j], UB[i-1][j] + T);
            if(j>0)   UB[i][j] = std::min(UB[i][j], UB[i][j-1] + T);
            if(i<m-1) UB[i][j] = std::min(UB[i][j], UB[i+1][j] - 1);
            if(j<h-1) UB[i][j] = std::min(UB[i][j], UB[i][j+1] - 1);
            if(oL!=LB[i][j] || oU!=UB[i][j]) changed=true;
        }
    }
}

// ------------- APT 多标签 DP（ages + 活跃列 + 自适应 K）-------------
CapTDAGResult ocCapTDAG::solve(int m,int h){
    LOGOPEN(cfg_.logfile, cfg_.log_append, cfg_.vlevel);

    const int n = m*h;
    int T = cfg_.T>0 ? cfg_.T : n;
    const uint64_t B = (uint64_t)h + 1u;

    // powB for base-(h+1)
    vector<uint64_t> powB(m);
    powB[m-1]=1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t t = (__uint128_t)powB[k+1]*(__uint128_t)B;
        if(t > (__uint128_t)ULLONG_MAX) throw std::runtime_error("state key overflow");
        powB[k] = (uint64_t)t;
    }

    LOG(1, [&](std::ofstream& o){
        o<<"[APT-ages+] m="<<m<<" h="<<h<<" n="<<n<<" T="<<T
         <<" v="<<cfg_.vlevel<<" K0="<<cfg_.max_labels_per_key<<"\n";
    });

    // windows
    vector<vector<int>> LB, UB;
    compute_windows_spanT(m,h,T,n,LB,UB);
    LOG(2, [&](std::ofstream& o){ o<<"[SpanT] window sanity check passed\n"; });

    // Label：相对龄值（0..T-1），未活跃维强制为 0
    struct Label{
        long long dist=0;
        uint64_t prev_key=0;
        int prev_row=-1;
        int prev_label_idx=-1;
        vector<uint16_t> r_age; // size=m（仅 a[i]>0 有意义）
        vector<uint16_t> c_age; // size=h（仅活跃列 J_active 有意义）
    };

    using LVec = vector<Label>;
    std::unordered_map<uint64_t, LVec> cur, nxt;

    auto build_active = [&](const vector<int>& a, vector<int>& J){
        J.clear();
        for(int i=1;i<m;++i) if(a[i]>0) J.push_back(a[i]);
        std::sort(J.begin(),J.end());
        J.erase(std::unique(J.begin(),J.end()), J.end());
    };

    const uint16_t AGE_CAP = (uint16_t)std::max(0, T-1); // 饱和到 T-1；T 表示下一步必超

    auto dominates = [&](const Label& B, const Label& A,
                         const vector<int>& a, const vector<int>& J)->bool{
        if(B.dist > A.dist) return false;
        bool strict = (B.dist < A.dist);
        for(int i=0;i<m;++i){
            if(a[i]==0) continue; // 未启动行不比较
            if(B.r_age[i] > A.r_age[i]) return false;
            if(B.r_age[i] < A.r_age[i]) strict = true;
        }
        for(int x : J){
            if(B.c_age[x] > A.c_age[x]) return false;
            if(B.c_age[x] < A.c_age[x]) strict = true;
        }
        return strict;
    };

    auto score = [&](const Label& L, const vector<int>& a, const vector<int>& J){
        long long s = 0;
        for(int i=0;i<m;++i) if(a[i]>0) s += L.r_age[i];
        for(int x : J) s += L.c_age[x];
        return std::pair<long long,long long>(L.dist, s); // 词典序：dist 主导，其次总龄值
    };

    int Kcur0 = std::max(1, cfg_.max_labels_per_key);
    const int KMAX = 4096;

    auto insert_label = [&](LVec& vec, Label&& L, int Kcur,
                            const vector<int>& a, const vector<int>& J,
                            size_t& pruned_dom, size_t& truncatedK){
        // 支配过滤
        for(const auto& e : vec) if(dominates(e, L, a, J)){ ++pruned_dom; return; }
        size_t w=0;
        for(size_t t=0;t<vec.size();++t){
            if(dominates(L, vec[t], a, J)){ ++pruned_dom; continue; }
            if(w!=t) vec[w] = std::move(vec[t]);
            ++w;
        }
        vec.resize(w);
        vec.push_back(std::move(L));

        // 截断：保留 Kcur 个最优 (dist, sum_active_ages)
        if((int)vec.size() > Kcur){
            std::nth_element(vec.begin(), vec.begin()+Kcur, vec.end(),
                [&](const Label& A, const Label& B){ return score(A,a,J) < score(B,a,J); });
            truncatedK += (vec.size() - Kcur);
            vec.resize(Kcur);
        }
    };

    // init
    uint64_t key0=0ull;
    {
        Label L0;
        L0.dist=0; L0.prev_key=0; L0.prev_row=-1; L0.prev_label_idx=-1;
        L0.r_age.assign(m, 0);
        L0.c_age.assign(h, 0);
        cur.emplace(key0, LVec{std::move(L0)});
    }

    // 保存每层的前沿（用于回溯）
    vector<std::unordered_map<uint64_t, LVec>> LV(n+1);

    for(int L=0; L<n; ++L){
        bool layer_ok=false;
        int retries=0;
        int Kcur = Kcur0;

        do{
            nxt.clear();
            nxt.reserve(std::max<size_t>(cur.size()*2, 16));

            size_t cand=0, blk_oc=0, blk_win=0, blk_tL=0, blk_tU=0, accepted=0, pruned_dom=0, truncatedK=0;

            for(auto &kv : cur){
                uint64_t key = kv.first;
                const LVec& labels = kv.second;

                vector<int> a(m);
                for(int i=0;i<m;++i) a[i] = digit_at(key,i,powB,B);

                vector<int> J; build_active(a, J); // 活跃列集合

                for(size_t li=0; li<labels.size(); ++li){
                    const Label& lab = labels[li];

                    for(int i=0;i<m;++i){
                        int j = a[i];
                        if(j>=h) continue;
                        ++cand;

                        // OC 非增：i>0 时需 j+1 ≤ a[i-1]
                        if(i>0 && j+1 > a[i-1]){ ++blk_oc; continue; }

                        int Lnext = L+1;

                        // 窗口必要：Lnext ∈ [LB,UB]
                        if(!(LB[i][j] <= Lnext && Lnext <= UB[i][j])){ ++blk_win; continue; }

                        // Δ≤T：ages 判定（age + 1 ≤ T ↔ age ≤ T-1）
                        if(j>0 && lab.r_age[i] >= (uint16_t)T){ ++blk_tL; continue; } // 左父
                        if(i>0 && lab.c_age[j] >= (uint16_t)T){ ++blk_tU; continue; } // 上父（j 属于 J_active）

                        uint64_t key2 = enc_inc(key, i, powB);
                        long long d1 = lab.dist + (long long)Lnext * weight_ij(i,j,m,h) * cfg_.dV;

                        // 更新 ages：仅活跃维 +1 饱和；本行/本列清零；未启动行置 0
                        vector<uint16_t> r2 = lab.r_age, c2 = lab.c_age;
                        for(int k=0;k<m;++k){
                            if(a[k]>0) r2[k] = sat_inc(r2[k], AGE_CAP);
                            else       r2[k] = 0;
                        }
                        for(int t=0;t<h;++t){
                            if(std::binary_search(J.begin(), J.end(), t)) c2[t] = sat_inc(c2[t], AGE_CAP);
                            else c2[t] = 0;
                        }
                        r2[i] = 0;
                        c2[j] = 0;
                        if(i==0){ // 新开一列 j=a0
                            int a0 = a[0];
                            if(a0 < h) c2[a0] = 0;
                        }

                        Label Lnew;
                        Lnew.dist=d1; Lnew.prev_key=key; Lnew.prev_row=i; Lnew.prev_label_idx=(int)li;
                        Lnew.r_age = std::move(r2);
                        Lnew.c_age = std::move(c2);

                        insert_label(nxt[key2], std::move(Lnew), Kcur, a, J, pruned_dom, truncatedK);
                        ++accepted;
                    }
                }
            }

            // 统计 + 自适应 K（不仅空前沿才扩容）
            size_t next_labels=0; for(auto &p: nxt) next_labels += p.second.size();
            double avg = nxt.empty()? 0.0 : (double)next_labels / (double)nxt.size();
            double pressure = accepted ? (double)truncatedK / (double)accepted : 0.0;

            if(cfg_.vlevel>=2){
                LOG(2, [&](std::ofstream& o){
                    o<<"[DP] L="<<L
                     <<" states="<<cur.size()<<" -> next="<<nxt.size()
                     <<" | Kcur="<<Kcur<<" retries="<<retries
                     <<" | pressure="<<pressure<<" avg/key="<<avg
                     <<" | cand="<<cand<<" truncK="<<truncatedK
                     <<" blk_oc="<<blk_oc<<" blk_win="<<blk_win
                     <<" blk_tL="<<blk_tL<<" blk_tU="<<blk_tU
                     <<"\n";
                });
            }

            // 触发规则：空前沿 或 截断压力大（≥0.3）或 avg/key 逼近 Kcur（≥0.9·Kcur）
            if( (nxt.empty() || pressure>=0.30 || (avg >= 0.90*Kcur)) && Kcur < KMAX ){
                int old = Kcur; Kcur = std::min(KMAX, Kcur*2); ++retries;
                LOG(1, [&](std::ofstream& o){ o<<"[DP] expand K "<<old<<" -> "<<Kcur<<" at L="<<L<<"\n"; });
                // 重做本层
            }else if(nxt.empty()){
                // 已到 KMAX 仍空前沿 ⇒ 确实无解
                throw std::runtime_error("No feasible next frontier (even after K growth)");
            }else{
                // 本层通过
                layer_ok = true;
                // 下一层的起始 K0 用当前 Kcur（避免层间振荡）
                Kcur0 = Kcur;
            }

        }while(!layer_ok);

        LV[L] = std::move(cur);
        cur.swap(nxt);
    }
    LV[n] = cur;

    // 终态挑最优标签
    uint64_t keyN = 0ull;
    for(int i=0;i<m;++i) keyN += powB[i]*(uint64_t)h;
    auto it = LV[n].find(keyN);
    if(it==LV[n].end() || it->second.empty()){
        throw std::runtime_error("No terminal state");
    }
    int best_idx=0;
    for(int t=1;t<(int)it->second.size();++t)
        if(it->second[t].dist < it->second[best_idx].dist) best_idx=t;

    // 回溯
    vector<vector<int>> y(m, vector<int>(h,0));
    uint64_t k = keyN; int li = best_idx;
    for(int L=n; L>=1; --L){
        const auto &mp = LV[L];
        auto hit = mp.find(k);
        if(hit==mp.end()) throw std::runtime_error("Reconstruct failed (key missing)");
        const vector<Label>& V = hit->second;
        if(li<0 || li>=(int)V.size()) li = 0; // 兜底

        const Label& curL = V[li];
        uint64_t pk = curL.prev_key;
        int ri = curL.prev_row;
        int aj = digit_at(pk, ri, powB, B);
        if(ri<0||ri>=m||aj<0||aj>=h) throw std::runtime_error("Reconstruct index OOR");

        y[ri][aj] = L;
        k = pk;
        li = curL.prev_label_idx;
        if(li<0) li=0;
    }

    long long total = hpwl_neighbors(y, cfg_.dV);
    int dmax = max_delta_adj(y);
    bool oc_ok = check_OC(y);

    LOG(1, [&](std::ofstream& o){
        o<<"[Verify] HPWL="<<total<<" maxΔ="<<dmax<<" OC="<<(oc_ok?"OK":"FAIL")<<"\n";
    });

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
