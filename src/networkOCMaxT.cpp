#include "networkOCMaxT.h"
#include <limits>
#include <fstream>
#include <stdexcept>
#include <algorithm>

using std::vector;
using std::unordered_map;
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
static inline int LB0(int i,int j){ return (i+1)*(j+1); }
}

// ======= 窗口闭包（与 lightOCMaxT 同构；仅为自包含重写） =======
void networkOCMaxT::compute_windows_spanT(
    int m,int h,int T, int n,
    vector<vector<int>>& LB,
    vector<vector<int>>& UB,
    int vlevel,
    const std::string& logpath,
    bool append)
{
    LOGOPEN2(logpath, append, vlevel);
    LB.assign(m, vector<int>(h, 0));
    UB.assign(m, vector<int>(h, 0));
    for(int i=0;i<m;++i)
        for(int j=0;j<h;++j){
            LB[i][j] = LB0(i,j);
            UB[i][j] = n - (m-1-i)*(h-1-j);
        }

    bool changed = true;
    int iter = 0, iter_max = m*h*6;
    while(changed && iter++ < iter_max){
        changed = false;
        // 父→子
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            int oldLB = LB[i][j], oldUB = UB[i][j];
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
            int oldLB = LB[i][j], oldUB = UB[i][j];
            if(i<m-1) LB[i][j] = std::max(LB[i][j], LB[i+1][j] - T);
            if(j<h-1) LB[i][j] = std::max(LB[i][j], LB[i][j+1] - T);
            if(i<m-1) UB[i][j] = std::min(UB[i][j], UB[i+1][j] - 1);
            if(j<h-1) UB[i][j] = std::min(UB[i][j], UB[i][j+1] - 1);
            if(LB[i][j]!=oldLB || UB[i][j]!=oldUB) changed = true;
        }
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            if(LB[i][j] > UB[i][j])
                throw std::runtime_error("Δ≤T infeasible after closure");
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
}

// ======= 主流程：完整状态 DAG 最短路 =======
OCMaxTResult networkOCMaxT::solve(int m, int h, int T_in){
    LOGOPEN2(cfg_.logfile, cfg_.log_append, cfg_.verbose_level);

    const int n_full = m*h;
    int T = T_in;
    if(T<=0) T=n_full;
    if(T>n_full) T=n_full;

    LOG2(1, [&](std::ofstream& o){
        o<<"[OC-SpanT-NET] n="<<n_full<<" T="<<T<<" v="<<cfg_.verbose_level
         <<" (full DAG shortest path)\n";
    });

    // 预计算 B 进制
    const uint64_t B = (uint64_t)h + 1u;
    vector<uint64_t> powB(m);
    powB[m-1] = 1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t tmp = (__uint128_t)powB[k+1]*(__uint128_t)B;
        if(tmp > (__uint128_t)std::numeric_limits<uint64_t>::max())
            throw std::runtime_error("state key overflow");
        powB[k]=(uint64_t)tmp;
    }
    auto get_level_from_key = [&](uint64_t key)->int{
        int s=0; for(int i=0;i<m;++i) s += digit_at(key,i,powB,B); return s;
    };

    // 窗口闭包
    vector<vector<int>> LB, UB;
    compute_windows_spanT(m,h,T,n_full,LB,UB,cfg_.verbose_level,cfg_.logfile,cfg_.log_append);
    LOG2(2, [&](std::ofstream& o){ o<<"[SpanT] window sanity check passed\n"; });

    struct Node { long long dist; uint64_t prev; int prev_i; };
    vector< unordered_map<uint64_t, Node> > layers(n_full+1);
    layers[0].reserve(1);
    layers[0].emplace(0ull, Node{0, 0ull, -1});
    uint64_t key_goal=0ull; for(int i=0;i<m;++i) key_goal += powB[i]*(uint64_t)h;

    // 逐层松弛（DAG 最短路）
    // …… 省略前置定义 ……
    for(int L=0; L<n_full; ++L){
    auto &cur = layers[L];
    auto &nxt = layers[L+1];
    nxt.reserve(std::max<size_t>(cur.size()*4, 16));

    size_t pruned_by_win=0, pruned_by_delta=0, gen_edges=0;

    for(const auto& kv : cur){
        uint64_t key = kv.first;
        long long d  = kv.second.dist;
        int a[64];
        for(int i=0;i<m;++i) a[i] = digit_at(key,i,powB,B);

        for(int i=0;i<m;++i){
            int ai = a[i];
            if(ai>=h) continue;
            if(i>0 && ai+1>a[i-1]) continue; // skyline 非增（OC）

            const int Lnext = L+1;

            // 窗口存在性
            if(!(LB[i][ai] <= Lnext && Lnext <= UB[i][ai])){ ++pruned_by_win; continue; }

            // === Δ≤T：用 min(UB_parent, L) 做存在性充分判定 ===
            if(ai>0){
                int parUB = std::min(UB[i][ai-1], L);
                if(Lnext > parUB + T){ ++pruned_by_delta; continue; }
            }
            if(i>0){
                int parUB = std::min(UB[i-1][ai], L);
                if(Lnext > parUB + T){ ++pruned_by_delta; continue; }
            }

            // 生成转移
            uint64_t key2 = encode_inc(key,i,powB);
            long long c = (long long)Lnext * weight_ij(i, ai, m, h) * cfg_.dV;
            long long nd = d + c;
            ++gen_edges;

            auto it = nxt.find(key2);
            if(it==nxt.end()){
                nxt.emplace(key2, Node{nd, key, i});
            }else if(nd < it->second.dist){
                it->second = Node{nd, key, i};
            }
        }
    }
    // ……（日志打印保持）……
}


    // 恢复路径
    auto itg = layers[n_full].find(key_goal);
    if(itg==layers[n_full].end()){
        LOG2(1, [&](std::ofstream& o){
            o<<"[NET] cannot reach goal state\n";
        });
        throw std::runtime_error("Goal unreachable");
    }
    vector<vector<int>> y(m, vector<int>(h,0));
    uint64_t k = key_goal;
    for(int L=n_full-1; L>=0; --L){
        const auto &nd = layers[L+1].at(k);
        int i = nd.prev_i;
        if(i<0) throw std::runtime_error("Reconstruction error");
        int ai_prev[64];
        for(int r=0;r<m;++r) ai_prev[r] = digit_at(nd.prev,r,powB,B);
        int j = ai_prev[i];
        y[i][j] = L+1;
        k = nd.prev;
        if(L==0) break;
    }

    // 评估
    long long HPWL=0;
    for(int i=0;i<m;++i)
        for(int j=0;j+1<h;++j) HPWL += llabs((long long)y[i][j+1]-y[i][j])*cfg_.dV;
    for(int i=0;i+1<m;++i)
        for(int j=0;j<h;++j) HPWL += llabs((long long)y[i+1][j]-y[i][j])*cfg_.dV;

    int dmax=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) dmax=std::max(dmax, y[i][j+1]-y[i][j]);
    for(int j=0;j<h;++j) for(int i=0;i+1<m;++i) dmax=std::max(dmax, y[i+1][j]-y[i][j]);
    if(dmax>T){
        LOG2(1, [&](std::ofstream& o){ o<<"[POST-NET] maxΔ="<<dmax<<" > T="<<T<<"\n"; });
        throw std::runtime_error("Post-check failed: max Δ > T");
    }
    LOG2(1, [&](std::ofstream& o){ o<<"[POST-NET] HPWL="<<HPWL<<" maxΔ="<<dmax<<" OK\n"; });

    OCMaxTResult R;
    R.m=m; R.h=h; R.n_full=n_full; R.T=T; R.dV=cfg_.dV;
    R.y_order=std::move(y);
    R.hpwl_full=HPWL;
    R.total_cost=HPWL;
    R.hpwl_prefix=0;
    R.oc_ok=true;
    return R;
}
