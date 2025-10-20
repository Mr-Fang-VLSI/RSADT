#include "networkOCMaxT.h"
#include "lightOCMaxT.h"    // 复用 OCMaxTResult、base 编码工具思想
#include <vector>
#include <unordered_map>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>

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

// --------- helpers ----------
static inline long long weight_ij(int i,int j,int m,int h){
    long long w=0;
    if(i==0)     w -= 1;
    if(i==m-1)   w += 1;
    if(j==0)     w -= 1;
    if(j==h-1)   w += 1;
    return w;
}
static inline int LB0(int i,int j){ return (i+1)*(j+1); }
static inline long long llabsll(long long x){ return x>=0?x:-x; }

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

// ====== 状态 Key 编码（a + R + C）======
// 为了简洁稳健，这里用字节串作为 key：8B(a_key) + 2B*m(R) + 2B*h(C)
static inline void append_u16(std::string& s, uint16_t v){
    char b[2]; b[0] = (char)(v & 0xFFu); b[1] = (char)((v>>8) & 0xFFu);
    s.append(b, b+2);
}
static inline uint16_t read_u16(const std::string& s, size_t off){
    return (uint16_t)((unsigned char)s[off] | ((unsigned char)s[off+1] << 8));
}
static inline void append_u64(std::string& s, uint64_t v){
    char b[8];
    for(int k=0;k<8;++k) b[k] = (char)((v>>(8*k)) & 0xFFull);
    s.append(b, b+8);
}
static inline uint64_t read_u64(const std::string& s, size_t off){
    uint64_t v=0;
    for(int k=7;k>=0;--k) v = (v<<8) | (unsigned char)s[off+k];
    return v;
}
static inline std::string make_key(uint64_t akey, const vector<int>& R, const vector<int>& C){
    std::string s; s.reserve(8 + 2*R.size() + 2*C.size());
    append_u64(s, akey);
    for(int x: R) append_u16(s, (uint16_t)x);
    for(int x: C) append_u16(s, (uint16_t)x);
    return s;
}
static inline void parse_key(const std::string& s, int m, int h,
                             uint64_t& akey, vector<int>& R, vector<int>& C){
    akey = read_u64(s, 0);
    R.resize(m); C.resize(h);
    size_t p=8;
    for(int i=0;i<m;++i){ R[i] = (int)read_u16(s, p); p+=2; }
    for(int j=0;j<h;++j){ C[j] = (int)read_u16(s, p); p+=2; }
}

// ====== 主流程：RCDC 完整状态 DAG 最短路 ======
OCMaxTResult networkOCMaxT::solve(int m, int h, int T_in){
    LOGOPEN2(cfg_.logfile, cfg_.log_append, cfg_.verbose_level);

    const int n_full = m*h;
    int T = T_in;
    if(T<=0) T=n_full;
    if(T>n_full) T=n_full;

    LOG2(1, [&](std::ofstream& o){
        o<<"[OC-SpanT-NET] n="<<n_full<<" T="<<T<<" v="<<cfg_.verbose_level
         <<" (RCDC full DAG shortest path)\n";
    });

    // base-(h+1) 编码
    const uint64_t B = (uint64_t)h + 1u;
    vector<uint64_t> powB(m);
    powB[m-1] = 1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t tmp = (__uint128_t)powB[k+1]*(__uint128_t)B;
        if(tmp > (__uint128_t)std::numeric_limits<uint64_t>::max())
            throw std::runtime_error("state key overflow");
        powB[k]=(uint64_t)tmp;
    }
    auto digit_at = [&](uint64_t key, int idx)->int{
        return (int)((key / powB[idx]) % B);
    };
    auto encode_inc = [&](uint64_t key, int idx)->uint64_t{
        return key + powB[idx];
    };

    // 窗口闭包
    vector<vector<int>> LB, UB;
    compute_windows_spanT(m,h,T,n_full,LB,UB,cfg_.verbose_level);

    // 初始状态
    uint64_t a0=0ull, aGoal=0ull;
    for(int i=0;i<m;++i) aGoal += powB[i]*(uint64_t)h;
    vector<int> R0(m, 0), C0(h, 0); // 0 表示“尚未定义”（未放置）
    std::string key0 = make_key(a0, R0, C0);

    struct Node { long long dist; std::string prev; int prev_i; };
    vector< unordered_map<std::string, Node> > layers(n_full+1);
    layers[0].reserve(1);
    layers[0].emplace(key0, Node{0, std::string(), -1});

    // 逐层最短路
    for(int L=0; L<n_full; ++L){
        auto &cur = layers[L];
        auto &nxt = layers[L+1];
        nxt.reserve(std::max<size_t>(cur.size()*4, 16));

        size_t pruned_win=0, pruned_delta=0, gen_edges=0;

        for(const auto& kv : cur){
            const std::string &K = kv.first;
            const Node &Nd = kv.second;

            // 解码
            uint64_t akey; vector<int> R(m), C(h);
            parse_key(K, m, h, akey, R, C);
            int a[64];
            for(int i=0;i<m;++i) a[i] = digit_at(akey, i);

            for(int i=0;i<m;++i){
                int j = a[i];
                if(j>=h) continue;
                if(i>0 && j+1>a[i-1]) continue; // skyline 非增（OC）

                int Lnext = L+1;

                // 窗口存在性（必要）
                if(!(LB[i][j] <= Lnext && Lnext <= UB[i][j])) { ++pruned_win; continue; }

                // === Δ≤T 精确检查（行/列父的“实际时间”） ===
                // 水平父：(i,j-1) 在同一行最后一次放置时间 R[i]
                if(j>0){
                    if(R[i]==0){ ++pruned_delta; continue; }        // 理论上不应发生；安全兜底
                    if(Lnext > R[i] + T){ ++pruned_delta; continue; }
                }
                // 竖直父：(i-1,j) 的实际时间就是当前列时钟 C[j]
                if(i>0){
                    if(C[j]==0){ ++pruned_delta; continue; }        // 竖直父尚不存在，非法
                    if(Lnext > C[j] + T){ ++pruned_delta; continue; }
                }

                // 生成后继状态
                uint64_t a2 = encode_inc(akey, i);
                vector<int> R2 = R, C2 = C;
                R2[i] = Lnext;
                C2[j] = Lnext;

                long long c = (long long)Lnext * weight_ij(i, j, m, h) * cfg_.dV;
                long long nd = Nd.dist + c;

                std::string K2 = make_key(a2, R2, C2);
                auto it = nxt.find(K2);
                if(it==nxt.end()){
                    nxt.emplace(std::move(K2), Node{nd, K, i});
                }else if(nd < it->second.dist){
                    it->second.dist = nd;
                    it->second.prev = K;
                    it->second.prev_i = i;
                }
                ++gen_edges;
            }
        }

        if(cfg_.progress && cfg_.verbose_level>=2 && ((L+1)%16==0)){
            LOG2(2, [&](std::ofstream& o){
                o<<"[NET-RCDC] level "<<L<<" states="<<cur.size()
                 <<" -> next="<<layers[L+1].size()
                 <<" | gen_edges="<<gen_edges
                 <<" pruned(win)="<<pruned_win
                 <<" pruned(Δ)="<<pruned_delta<<"\n";
            });
        }
        if(layers[L+1].empty()){
            LOG2(1, [&](std::ofstream& o){
                o<<"[NET-RCDC] dead at level "<<(L+1)<<", no feasible transitions\n";
            });
            throw std::runtime_error("No feasible state at some level (RCDC)");
        }
    }

    // 在最后一层选择 a=Goal 的最优状态
    long long best = std::numeric_limits<long long>::max();
    std::string Kbest;
    for(const auto& kv : layers[n_full]){
        uint64_t akey; vector<int> R(m), C(h);
        parse_key(kv.first, m, h, akey, R, C);
        if(akey!=aGoal) continue;
        if(kv.second.dist < best){ best = kv.second.dist; Kbest = kv.first; }
    }
    if(Kbest.empty()){
        LOG2(1, [&](std::ofstream& o){ o<<"[NET-RCDC] cannot reach goal state\n"; });
        throw std::runtime_error("Goal unreachable (RCDC)");
    }

    // 复原 y
    vector<vector<int>> y(m, vector<int>(h, 0));
    std::string K = Kbest;
    for(int L=n_full-1; L>=0; --L){
        const Node &Nd = layers[L+1].at(K);
        // 找出是哪一行增加
        uint64_t akey; vector<int> R(m), C(h);
        parse_key(K, m, h, akey, R, C);
        // 前驱 a
        uint64_t akey_prev; vector<int> Rprev(m), Cprev(h);
        parse_key(Nd.prev, m, h, akey_prev, Rprev, Cprev);
        int sel_i=-1, sel_j=-1;
        for(int i=0;i<m;++i){
            int ai_prev = (int)((akey_prev / powB[i]) % B);
            int ai_now  = (int)((akey      / powB[i]) % B);
            if(ai_now == ai_prev+1){ sel_i = i; sel_j = ai_prev; break; }
        }
        if(sel_i<0) throw std::runtime_error("Reconstruct mismatch");
        y[sel_i][sel_j] = L+1;
        K = Nd.prev;
        if(L==0) break;
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
        throw std::runtime_error("Post-check failed: max Δ > T (RCDC)"); // 理论不应触发
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
