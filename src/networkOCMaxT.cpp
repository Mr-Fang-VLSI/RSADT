#include "networkOCMaxT.h"
#include "lightOCMaxT.h"

#include <vector>
#include <unordered_map>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdint>

using std::vector;
using std::unordered_map;

// ---------------- Logger ----------------
namespace {
std::ofstream g_log2;
int g_v2 = 1;
inline void LOGOPEN2(const std::string& path, bool append, int v) {
    if (g_log2.is_open()) return;
    g_log2.open(path, append ? (std::ios::out|std::ios::app) : (std::ios::out|std::ios::trunc));
    g_v2 = v;
}
template<typename F>
inline void LOG2(int lvl, F&& f){
    if (g_v2 >= lvl && g_log2.is_open()) { f(g_log2); g_log2.flush(); }
}
}

// ---------- helpers ----------
static inline long long weight_ij(int i,int j,int m,int h){
    long long w=0;
    if(i==0)     w -= 1;
    if(i==m-1)   w += 1;
    if(j==0)     w -= 1;
    if(j==h-1)   w += 1;
    return w;
}
static inline long long llabsll(long long x){ return x>=0?x:-x; }
static inline int LB0(int i,int j){ return (i+1)*(j+1); }

static inline uint16_t sat_inc(uint16_t v, uint16_t T){ return (v<T)?(uint16_t)(v+1):T; }

// -------- Δ≤T: 窗口闭包（存在性必要 + OC 单调性）---------
static void compute_windows_spanT(
    int m,int h,int T, int n,
    vector<vector<int>>& LB,
    vector<vector<int>>& UB,
    int vlevel)
{
    LB.assign(m, vector<int>(h, 0));
    UB.assign(m, vector<int>(h, 0));
    for(int i=0;i<m;++i)
        for(int j=0;j<h;++j){
            LB[i][j] = LB0(i,j);                           // 基本最早（OC）
            UB[i][j] = n - (m-1-i)*(h-1-j);               // 基本最晚（OC）
        }

    bool changed = true;
    int iter = 0, iter_max = m*h*6;
    while(changed && iter++ < iter_max){
        changed = false;
        // 父→子
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            int oldLB=LB[i][j], oldUB=UB[i][j];
            if(i>0)   LB[i][j] = std::max(LB[i][j], LB[i-1][j] + 1);
            if(j>0)   LB[i][j] = std::max(LB[i][j], LB[i][j-1] + 1);
            if(i>0)   UB[i][j] = std::min(UB[i][j], UB[i-1][j] + T);
            if(j>0)   UB[i][j] = std::min(UB[i][j], UB[i][j-1] + T);
            if(i<m-1) UB[i][j] = std::min(UB[i][j], UB[i+1][j] - 1);
            if(j<h-1) UB[i][j] = std::min(UB[i][j], UB[i][j+1] - 1);
            if(LB[i][j]!=oldLB || UB[i][j]!=oldUB) changed = true;
        }
        // 子→父
        for(int i=m-1;i>=0;--i) for(int j=h-1;j>=0;--j){
            int oldLB=LB[i][j], oldUB=UB[i][j];
            if(i<m-1) LB[i][j] = std::max(LB[i][j], LB[i+1][j] - T);
            if(j<h-1) LB[i][j] = std::max(LB[i][j], LB[i][j+1] - T);
            if(i<m-1) UB[i][j] = std::min(UB[i][j], UB[i+1][j] - 1);
            if(j<h-1) UB[i][j] = std::min(UB[i][j], UB[i][j+1] - 1);
            if(LB[i][j]!=oldLB || UB[i][j]!=oldUB) changed = true;
        }
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            if(LB[i][j] > UB[i][j]) throw std::runtime_error("Δ≤T infeasible after closure");
        }
    }
    if (vlevel>=3){
        int LBmin=1e9, LBmax=-1e9, UBmin=1e9, UBmax=-1e9;
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            LBmin = std::min(LBmin, LB[i][j]);
            LBmax = std::max(LBmax, LB[i][j]);
            UBmin = std::min(UBmin, UB[i][j]);
            UBmax = std::max(UBmax, UB[i][j]);
        }
        LOG2(3, [&](std::ofstream& o){
            o<<"[SpanT] LB range=["<<LBmin<<","<<LBmax<<"], UB range=["
              <<UBmin<<","<<UBmax<<"]\n";
        });
    }
    LOG2(2, [&](std::ofstream& o){ o<<"[SpanT] window sanity check passed\n"; });
}

// ------------- base-(h+1) 编解码 -------------
static void build_powB(int m, int h, vector<uint64_t>& powB, uint64_t& B){
    B = (uint64_t)h + 1u;
    powB.resize(m);
    powB[m-1] = 1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t tmp = (__uint128_t)powB[k+1]*(__uint128_t)B;
        if(tmp > (__uint128_t)std::numeric_limits<uint64_t>::max())
            throw std::runtime_error("state key overflow");
        powB[k]=(uint64_t)tmp;
    }
}
static inline int digit_at(uint64_t key, int idx, const vector<uint64_t>& powB, uint64_t B){
    return (int)((key / powB[idx]) % B);
}
static inline uint64_t encode_inc(uint64_t key, int idx, const vector<uint64_t>& powB){
    return key + powB[idx];
}

// ----------- 按天际线分桶 + 主导裁剪 ------------
struct NDState {
    long long dist{0};
    vector<uint16_t> r_age; // size = m；未启动行规范化为 0
    vector<uint16_t> c_age; // size = h；只在前缀 [0,a[0]) 上有效
    // 回溯
    uint64_t prev_akey{0};
    int      prev_bucket{-1};
    int      prev_state{-1};
    int      prev_row{-1};
};
struct Bucket {
    uint64_t akey{0};
    vector<NDState> states; // 同一 akey 的非支配集
};
struct Layer {
    vector<Bucket> buckets;
    unordered_map<uint64_t,int> a2idx; // akey -> buckets 下标
};

// 主导：同 a 下，行只比较 a[i]>0 的行；列只比较前缀 j<a0
static inline bool dominates(const NDState& ex, const NDState& cand,
                             const vector<int>& a, int a0)
{
    if(ex.dist > cand.dist) return false; // 代价必须不大于
    const int m = (int)a.size();
    for(int i=0;i<m;++i){
        if(a[i]==0) continue;                 // 未启动行不比较
        if(ex.r_age[i] > cand.r_age[i]) return false;
    }
    for(int j=0;j<a0;++j){                    // 只比较前缀
        if(ex.c_age[j] > cand.c_age[j]) return false;
    }
    return true;
}

static bool insert_nondominated(vector<NDState>& vec, NDState cand,
                                const vector<int>& a, int a0,
                                int cap_per_bucket, bool enable_dominance)
{
    if(enable_dominance){
        for(const auto& ex : vec){
            if(dominates(ex, cand, a, a0)) return false; // cand 被支配
        }
        int w=0;
        for(int k=0;k<(int)vec.size();++k){
            if(dominates(cand, vec[k], a, a0)) continue; // 丢弃被支配者
            vec[w++] = std::move(vec[k]);
        }
        vec.resize(w);
    }
    vec.push_back(std::move(cand));

    // 可选：前沿上限（稳定内存；0=不开）
    if(cap_per_bucket>0 && (int)vec.size()>cap_per_bucket){
        std::nth_element(vec.begin(), vec.begin()+cap_per_bucket, vec.end(),
                         [](const NDState& A, const NDState& B){ return A.dist < B.dist; });
        vec.resize(cap_per_bucket);
    }
    return true;
}

// ====== 主流程：RCDC + 相对龄值（固定长度 c_age） + 主导裁剪 ======
OCMaxTResult networkOCMaxT::solve(int m, int h, int T_in){
    LOGOPEN2(cfg_.logfile, cfg_.log_append, cfg_.verbose_level);

    const int n_full = m*h;
    int T = T_in;
    if(T<=0) T=n_full;
    if(T>n_full) T=n_full;

    LOG2(1, [&](std::ofstream& o){
        o<<"[OC-SpanT-NET] n="<<n_full<<" T="<<T<<" v="<<cfg_.verbose_level
         <<" (RCDC ages + dominance, fixed c_age)\n";
    });

    // base 编码
    vector<uint64_t> powB; uint64_t B=0;
    build_powB(m,h,powB,B);
    auto decode_a = [&](uint64_t akey)->vector<int>{
        vector<int> a(m);
        for(int i=0;i<m;++i) a[i] = digit_at(akey,i,powB,B);
        return a;
    };

    // 窗口闭包
    vector<vector<int>> LB, UB;
    compute_windows_spanT(m,h,T,n_full,LB,UB,cfg_.verbose_level);

    // 起点
    uint64_t a0key=0ull, aGoal=0ull;
    for(int i=0;i<m;++i) aGoal += powB[i]*(uint64_t)h;

    vector<Layer> layers; layers.reserve(n_full+1);
    {
        Layer L0; L0.buckets.reserve(1); L0.a2idx.reserve(1);
        Bucket b0; b0.akey=a0key;
        NDState s0;
        s0.dist=0;
        s0.r_age.assign(m, 0); // 未启动行置 0
        s0.c_age.assign(h, 0); // 固定长度；前缀长度 a[0]=0，此时都视为未启用
        b0.states.push_back(std::move(s0));
        L0.buckets.push_back(std::move(b0));
        L0.a2idx.emplace(a0key,0);
        layers.push_back(std::move(L0));
    }

    // 逐层扩展
    for(int L=0; L<n_full; ++L){
        const Layer& cur = layers.back();
        Layer nxt; nxt.buckets.reserve(cur.buckets.size()*2+16);
        nxt.a2idx.reserve(cur.buckets.size()*2+16);

        size_t states_cur=0, states_nxt=0, gen_edges=0;
        size_t pruned_win=0, pruned_delta=0, inserted=0;

        for(int b=0; b<(int)cur.buckets.size(); ++b){
            const Bucket& BK = cur.buckets[b];
            vector<int> a = decode_a(BK.akey);
            const int a0 = a[0];

            for(int s=0; s<(int)BK.states.size(); ++s){
                const NDState& S = BK.states[s];
                ++states_cur;

                // 每一行尝试放置
                for(int i=0;i<m;++i){
                    int j = a[i];
                    if(j>=h) continue;
                    if(i>0 && j+1 > a[i-1]) continue; // skyline 非增（OC）

                    int Lnext = L+1;

                    // 窗口存在性（必要）
                    if(!(LB[i][j] <= Lnext && Lnext <= UB[i][j])){ ++pruned_win; continue; }

                    // Δ≤T（相对龄值）：前缀内判断
                    if(j>0 && S.r_age[i] >= (uint16_t)T){ ++pruned_delta; continue; }
                    if(i>0){
                        if(j >= a0){ ++pruned_delta; continue; }               // 尚未启用该列
                        if(S.c_age[j] >= (uint16_t)T){ ++pruned_delta; continue; }
                    }

                    // 生成后继
                    vector<int> a2 = a; a2[i] = j+1;
                    uint64_t a2key = encode_inc(BK.akey, i, powB);
                    const int a0_new = a2[0];

                    // 更新龄值（相对、截到 T）
                    vector<uint16_t> r2 = S.r_age;
                    for(int k=0;k<m;++k){
                        r2[k] = sat_inc(r2[k], (uint16_t)T);
                    }
                    r2[i] = 0;                            // 本行重置
                    for(int k=0;k<m;++k) if(a2[k]==0) r2[k]=0;  // 未启动行归零（规范化）

                    vector<uint16_t> c2 = S.c_age;        // 固定长度 h
                    // 仅对已启用前缀 [0,a0) 增龄
                    for(int t=0;t<a0;++t) c2[t] = sat_inc(c2[t], (uint16_t)T);
                    if(i==0){
                        // 新增列（位置 a0）
                        if(a0_new==a0+1 && a0 < (int)c2.size()) c2[a0] = 0;
                    }else{
                        // 竖直父列 j 必在前缀内
                        c2[j] = 0;
                    }

                    long long edge_cost = (long long)Lnext * weight_ij(i, j, m, h) * cfg_.dV;
                    NDState NS;
                    NS.dist = S.dist + edge_cost;
                    NS.r_age = std::move(r2);
                    NS.c_age = std::move(c2);
                    NS.prev_akey = BK.akey;
                    NS.prev_bucket = b;
                    NS.prev_state  = s;
                    NS.prev_row    = i;

                    // 放入下一层的对应 bucket（按 a2key）
                    auto it = nxt.a2idx.find(a2key);
                    int bx = -1;
                    if(it==nxt.a2idx.end()){
                        bx = (int)nxt.buckets.size();
                        Bucket nb; nb.akey = a2key;
                        nxt.buckets.push_back(std::move(nb));
                        nxt.a2idx.emplace(a2key, bx);
                    }else bx = it->second;

                    if(insert_nondominated(nxt.buckets[bx].states,
                                           std::move(NS), a2, a0_new,
                                           cfg_.cap_per_bucket,
                                           cfg_.enable_dominance))
                    {
                        ++inserted;
                    }
                    ++gen_edges;
                } // for i
            } // for s
        } // for buckets

        for(const auto& b : nxt.buckets) states_nxt += b.states.size();

        if(cfg_.progress && ((cfg_.verbose_level>=3) ||
                             (cfg_.verbose_level>=2 && ((L+1)%std::max(1,cfg_.log_every_levels)==0)))){
            LOG2(2, [&](std::ofstream& o){
                o<<"[NET-RCDC] level "<<L
                 <<" buckets="<<cur.buckets.size()<<" -> "<<nxt.buckets.size()
                 <<" | states="<<states_cur<<" -> "<<states_nxt
                 <<" | gen_edges="<<gen_edges
                 <<" pruned(win)="<<pruned_win
                 <<" pruned(Δ)="<<pruned_delta
                 <<" inserted="<<inserted<<"\n";
            });
        }
        if(nxt.buckets.empty()){
            LOG2(1, [&](std::ofstream& o){
                o<<"[NET-RCDC] dead at level "<<(L+1)<<", no feasible transitions\n";
            });
            throw std::runtime_error("No feasible state at some level (RCDC-ages/fixed)");
        }
        layers.push_back(std::move(nxt));
    } // for L

    // 选择 a=Goal 的最优终态
    const Layer& Last = layers.back();
    auto itg = Last.a2idx.find( (uint64_t)aGoal );
    if(itg==Last.a2idx.end() || Last.buckets[itg->second].states.empty()){
        LOG2(1, [&](std::ofstream& o){ o<<"[NET-RCDC] cannot reach goal state\n"; });
        throw std::runtime_error("Goal unreachable (RCDC-ages/fixed)");
    }
    int best_idx = 0;
    const vector<NDState>& finals = Last.buckets[itg->second].states;
    for(int k=1;k<(int)finals.size();++k){
        if(finals[k].dist < finals[best_idx].dist) best_idx = k;
    }

    // 回溯 y
    vector<vector<int>> y(m, vector<int>(h,0));
    int Lcur = n_full;
    uint64_t cur_akey = aGoal;
    int cur_bucket = itg->second;
    int cur_state  = best_idx;

    auto decode_digit = [&](uint64_t akey, int idx)->int{
        return digit_at(akey, idx, powB, B);
    };

    while(Lcur>0){
        const NDState& S = layers[Lcur].buckets[cur_bucket].states[cur_state];
        int i = S.prev_row;
        if(i<0) throw std::runtime_error("Reconstruct error: invalid prev");

        const uint64_t prev_akey = S.prev_akey;
        int j = decode_digit(prev_akey, i); // 放置前该行的下一列索引
        y[i][j] = Lcur;

        cur_akey  = prev_akey;
        cur_bucket= S.prev_bucket;
        cur_state = S.prev_state;
        --Lcur;
    }

    // 评估与校验
    long long HPWL=0;
    for(int i=0;i<m;++i)
        for(int j=0;j+1<h;++j) HPWL += llabsll((long long)y[i][j+1]-y[i][j])*cfg_.dV;
    for(int i=0;i+1<m;++i)
        for(int j=0;j<h;++j)   HPWL += llabsll((long long)y[i+1][j]-y[i][j])*cfg_.dV;

    int dmax=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) dmax=std::max(dmax, y[i][j+1]-y[i][j]);
    for(int j=0;j<h;++j) for(int i=0;i+1<m;++i) dmax=std::max(dmax, y[i+1][j]-y[i][j]);

    if(dmax>T){
        LOG2(1, [&](std::ofstream& o){ o<<"[POST-RCDC] maxΔ="<<dmax<<" > T="<<T<<"\n"; });
        throw std::runtime_error("Post-check failed: max Δ > T (RCDC-ages/fixed)");
    }
    LOG2(1, [&](std::ofstream& o){ o<<"[POST-RCDC] HPWL="<<HPWL<<" maxΔ="<<dmax<<" OK\n"; });

    OCMaxTResult Rret;
    Rret.m=m; Rret.h=h; Rret.n_full=n_full; Rret.T=T; Rret.dV=cfg_.dV;
    Rret.y_order=std::move(y);
    Rret.hpwl_full=HPWL;
    Rret.total_cost=HPWL;
    Rret.hpwl_prefix=0;
    Rret.oc_ok=true;
    return Rret;
}
