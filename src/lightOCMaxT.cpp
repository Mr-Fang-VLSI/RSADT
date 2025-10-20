#include "lightOCMaxT.h"
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <fstream>
#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;
using std::unordered_map;

// ---------------- Logger ----------------
namespace {
std::ofstream g_log;
int g_vl = 1;
inline void LOGOPEN(const std::string& path, bool append, int v) {
    if (g_log.is_open()) return;
    g_log.open(path, append ? (std::ios::out|std::ios::app) : (std::ios::out|std::ios::trunc));
    g_vl = v;
}
template<typename F>
inline void LOG(int lvl, F&& f){
    if (g_vl >= lvl && g_log.is_open()) { f(g_log); g_log.flush(); }
}
}

// -------------- helpers -----------------
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

// 评估
long long lightOCMaxT::hpwl_full_neighbors(const vector<vector<int>>& y, long long dV){
    const int m = (int)y.size(), h = (int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i)
        for(int j=0;j+1<h;++j) S += llabsll((long long)y[i][j+1]-y[i][j])*dV;
    for(int i=0;i+1<m;++i)
        for(int j=0;j<h;++j)   S += llabsll((long long)y[i+1][j]-y[i][j])*dV;
    return S;
}
long long lightOCMaxT::hpwl_prefix_neighbors(const vector<vector<int>>& y, long long dV){
    const int m = (int)y.size(), h = (int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i)
        for(int j=0;j+1<h;++j){
            int a=y[i][j], b=y[i][j+1];
            if(a>0 && b>0) S += llabsll((long long)b-a)*dV;
        }
    for(int i=0;i+1<m;++i)
        for(int j=0;j<h;++j){
            int a=y[i][j], b=y[i+1][j];
            if(a>0 && b>0) S += llabsll((long long)b-a)*dV;
        }
    return S;
}
bool lightOCMaxT::check_OC_prefix(const vector<vector<int>>& y){
    const int m = (int)y.size(), h=(int)y[0].size();
    for(int i1=0;i1<m;++i1)
        for(int j1=0;j1<h;++j1){
            int v1 = y[i1][j1];
            if(v1<=0) continue;
            for(int i2=i1;i2<m;++i2)
                for(int j2=j1;j2<h;++j2){
                    int v2 = y[i2][j2];
                    if(v2<=0) continue;
                    if(v1 > v2) return false;
                }
        }
    return true;
}

struct KV { uint64_t key; long long dist; };

// ===== Δ≤T: 差分约束闭包（存在性必要 + OC 单调性）=====
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

        // ---------- 父→子：存在性必要 + OC ----------
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            int oldLB = LB[i][j], oldUB = UB[i][j];

            // OC：子 ≥ 父+1
            if(i>0)   LB[i][j] = std::max(LB[i][j], LB[i-1][j] + 1);
            if(j>0)   LB[i][j] = std::max(LB[i][j], LB[i][j-1] + 1);

            // Δ≤T（存在性必要）：子 ≤ 父的最大上界 + T
            if(i>0)   UB[i][j] = std::min(UB[i][j], UB[i-1][j] + T);
            if(j>0)   UB[i][j] = std::min(UB[i][j], UB[i][j-1] + T);

            // OC 单调性（更紧致）
            if(i<m-1) UB[i][j] = std::min(UB[i][j], UB[i+1][j] - 1);
            if(j<h-1) UB[i][j] = std::min(UB[i][j], UB[i][j+1] - 1);

            if(LB[i][j]!=oldLB || UB[i][j]!=oldUB) changed = true;
        }

        // ---------- 子→父：存在性必要 + OC ----------
        for(int i=m-1;i>=0;--i) for(int j=h-1;j>=0;--j){
            int oldLB = LB[i][j], oldUB = UB[i][j];

            // Δ≤T 反向：父 ≥ 子 - T（存在性必要）
            if(i<m-1) LB[i][j] = std::max(LB[i][j], LB[i+1][j] - T);
            if(j<h-1) LB[i][j] = std::max(LB[i][j], LB[i][j+1] - T);

            // OC 反向：父 ≤ 子 - 1
            if(i<m-1) UB[i][j] = std::min(UB[i][j], UB[i+1][j] - 1);
            if(j<h-1) UB[i][j] = std::min(UB[i][j], UB[i][j+1] - 1);

            if(LB[i][j]!=oldLB || UB[i][j]!=oldUB) changed = true;
        }

        // 一致性检查
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            if(LB[i][j] > UB[i][j]) {
                LOG(1, [&](std::ofstream& o){
                    o<<"[SpanT] infeasible at ("<<i<<","<<j<<") LB="<<LB[i][j]
                     <<" UB="<<UB[i][j]<<"\n";
                });
                throw std::runtime_error("Δ≤T infeasible after closure");
            }
        }
    }

#ifndef NDEBUG
    // Sanity（对应采用的规则）
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        if(i>0){
            assert(LB[i][j] >= LB[i-1][j] + 1);  // OC
            assert(UB[i][j] <= UB[i-1][j] + T);  // Δ≤T 存在性
        }
        if(j>0){
            assert(LB[i][j] >= LB[i][j-1] + 1);
            assert(UB[i][j] <= UB[i][j-1] + T);
        }
        assert(LB[i][j] <= UB[i][j]);
    }
    LOG(2, [&](std::ofstream& o){ o<<"[SpanT] window sanity check passed\n"; });
#else
    LOG(2, [&](std::ofstream& o){ o<<"[SpanT] window closure finished\n"; });
#endif

    if (vlevel>=3){
        int LBmin=1e9, LBmax=-1e9, UBmin=1e9, UBmax=-1e9;
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            LBmin = std::min(LBmin, LB[i][j]);
            LBmax = std::max(LBmax, LB[i][j]);
            UBmin = std::min(UBmin, UB[i][j]);
            UBmax = std::max(UBmax, UB[i][j]);
        }
        LOG(3, [&](std::ofstream& o){
            o<<"[SpanT] LB range=["<<LBmin<<","<<LBmax<<"], UB range=["
              <<UBmin<<","<<UBmax<<"]\n";
        });
    }
}

// --------- cut 兼容性检查：mid 状态能否保证所有跨切分边存在满足 1..T 的差分 ----------
static bool cut_state_T_compatible(
    int m,int h,int lmid,const vector<int>& a,int T,
    const vector<vector<int>>& LB, const vector<vector<int>>& UB)
{
    auto ok_pair = [&](int pi,int pj,int ci,int cj)->bool{
        int Pmin = LB[pi][pj];
        int Pmax = std::min(UB[pi][pj], lmid);
        int Cmin = std::max(LB[ci][cj], lmid+1);
        int Cmax = UB[ci][cj];
        if(Pmin>Pmax || Cmin>Cmax) return false;
        // 存在 r_p ∈ [Pmin,Pmax], r_c ∈ [Cmin,Cmax] 使 1 ≤ r_c - r_p ≤ T
        return (Cmax >= Pmin + 1) && (Cmin <= Pmax + T);
    };

    // 水平跨边
    for(int i=0;i<m;++i){
        int ai = a[i];
        if(0<ai && ai<h){
            if(!ok_pair(i, ai-1, i, ai)) return false;
        }
    }
    // 竖直跨边：当 a[i] > a[i+1]，j in [a[i+1], a[i)-1]
    for(int i=0;i+1<m;++i){
        if(a[i] > a[i+1]){
            for(int j=a[i+1]; j<=a[i]-1; ++j){
                if(!ok_pair(i, j, i+1, j)) return false;
            }
        }
    }
    return true;
}

// -------- 两层 DP（正向），加入 Δ≤T 的动态充分约束（基于父 LB） --------
static void forward_to_level_map(
    int m,int h, uint64_t key0,int l0,int steps, long long dV, int Tcap,
    const vector<vector<int>>* LBwin,
    const vector<vector<int>>* UBwin,
    const vector<uint64_t>& powB, uint64_t B,
    unordered_map<uint64_t,long long>& out,
    bool progress, int omp_threads, int vlevel)
{
    vector<KV> cur; cur.reserve(4096);
    cur.push_back({key0, 0});
    for(int lvl=0; lvl<steps; ++lvl){
        const int Lnext = l0 + lvl + 1;
        size_t est = std::max<size_t>(cur.size()*4, 16);
        vector<KV> next; next.reserve(est);
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(est);

#ifdef _OPENMP
        int TH = omp_threads>0? omp_threads : omp_get_max_threads();
        vector< vector<KV> > nxt_local(TH);
        vector< unordered_map<uint64_t,int> > id_local(TH);
        for(int t=0;t<TH;++t){ nxt_local[t].reserve(est/TH+16); id_local[t].reserve(est/TH+16); }

#pragma omp parallel num_threads(TH)
        {
            int tid=0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            auto &nl = nxt_local[tid];
            auto &il = id_local[tid];

#pragma omp for schedule(static)
            for(int u=0; u<(int)cur.size(); ++u){
                const uint64_t ku = cur[u].key;
                const long long du= cur[u].dist;
                int ai[32]; for(int i=0;i<m;++i) ai[i] = lightOCMaxT::digit_at(ku,i,powB,B);

                for(int i=0;i<m;++i){
                    int aii=ai[i];
                    if(aii>=h) continue;
                    if(i>0 && aii+1>ai[i-1]) continue;  // 非增约束（OC）

                    // 静态窗口存在性检查
                    if(LBwin && UBwin){
                        int lb = (*LBwin)[i][aii];
                        int ub = (*UBwin)[i][aii];
                        if(!(lb <= Lnext && Lnext <= ub)) continue;
                    }
                    // Δ≤T 动态充分约束（用父的 LB）
                    // Δ≤T 动态存在性充分：用 min(UB_parent, L) 作为父能取到的最晚层
if(LBwin && UBwin){
    if(aii>0){
        int parUB = std::min((*UBwin)[i][aii-1], Lnext-1); // Lnext-1 就是当前 L
        if(Lnext > parUB + Tcap) continue; // 左父
    }
    if(i>0){
        int parUB = std::min((*UBwin)[i-1][aii], Lnext-1);
        if(Lnext > parUB + Tcap) continue; // 上父
    }
}


                    uint64_t keyNew = lightOCMaxT::encode_digit_inc(ku,i,powB);
                    long long c = (long long)(Lnext) * weight_ij(i, aii, m, h) * dV;
                    long long nd = du + c;

                    auto it = il.find(keyNew);
                    if(it==il.end()){
                        int v=(int)nl.size();
                        nl.push_back({keyNew, nd});
                        il.emplace(keyNew, v);
                    }else{
                        int v = it->second;
                        if(nd<nl[v].dist) nl[v].dist=nd;
                    }
                }
            }
        }
        for(size_t t=0;t<nxt_local.size();++t){
            for(auto &kvpair : nxt_local[t]){
                auto it = nxt_id.find(kvpair.key);
                if(it==nxt_id.end()){
                    int v=(int)next.size();
                    next.push_back(kvpair);
                    nxt_id.emplace(kvpair.key, v);
                }else{
                    int v = it->second;
                    if(kvpair.dist<next[v].dist) next[v].dist=kvpair.dist;
                }
            }
            vector<KV>().swap(nxt_local[t]); unordered_map<uint64_t,int>().swap(id_local[t]);
        }
#else
        for(int u=0; u<(int)cur.size(); ++u){
            const uint64_t ku = cur[u].key;
            const long long du= cur[u].dist;
            int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCMaxT::digit_at(ku,i,powB,B);

            for(int i=0;i<m;++i){
                int aii=ai[i];
                if(aii>=h) continue;
                if(i>0 && aii+1>ai[i-1]) continue;

                if(LBwin && UBwin){
                    int lb = (*LBwin)[i][aii];
                    int ub = (*UBwin)[i][aii];
                    if(!(lb <= Lnext && Lnext <= ub)) continue;
                }
                if(LBwin){
                    if(aii>0 && Lnext > (*LBwin)[i][aii-1] + Tcap) continue;
                    if(i>0  && Lnext > (*LBwin)[i-1][aii] + Tcap) continue;
                }

                uint64_t keyNew = lightOCMaxT::encode_digit_inc(ku,i,powB);
                long long c = (long long)(Lnext) * weight_ij(i, aii, m, h) * dV;
                long long nd = du + c;

                auto it = nxt_id.find(keyNew);
                if(it==nxt_id.end()){
                    int v=(int)next.size();
                    next.push_back({keyNew, nd});
                    nxt_id.emplace(keyNew, v);
                }else{
                    int v = it->second;
                    if(nd<next[v].dist) next[v].dist=nd;
                }
            }
        }
#endif
        if(progress && (vlevel>=2) && (Lnext%16==0)){
            LOG(2, [&](std::ofstream& o){
                o<<"[DP-fwd] level "<<Lnext-1<<" states="<<cur.size()<<" -> next="<<next.size()<<"\n";
            });
        }
        cur.swap(next);
    }
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv : cur) out.emplace(kv.key, kv.dist);
}

// -------- 两层 DP（反向），不加动态 Δ 充分约束（避免过紧），仅做静态窗口 + OC ----------
static void backward_to_level_map(
    int m,int h, uint64_t keyN,int lN,int steps, long long dV, int /*Tcap*/,
    const vector<vector<int>>* LBwin,
    const vector<vector<int>>* UBwin,
    const vector<uint64_t>& powB, uint64_t B,
    unordered_map<uint64_t,long long>& out,
    bool progress, int omp_threads, int vlevel)
{
    vector<KV> cur; cur.reserve(4096);
    cur.push_back({keyN, 0});
    for(int t=0; t<steps; ++t){
        const int Lcur = lN - t;   // 正向 rank
        size_t est = std::max<size_t>(cur.size()*4, 16);
        vector<KV> next; next.reserve(est);
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(est);

#ifdef _OPENMP
        int TH = omp_threads>0? omp_threads : omp_get_max_threads();
        vector< vector<KV> > nxt_local(TH);
        vector< unordered_map<uint64_t,int> > id_local(TH);
        for(int k=0;k<TH;++k){ nxt_local[k].reserve(est/TH+16); id_local[k].reserve(est/TH+16); }

#pragma omp parallel num_threads(TH)
        {
            int tid=0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            auto &nl = nxt_local[tid];
            auto &il = id_local[tid];
#pragma omp for schedule(static)
            for(int u=0; u<(int)cur.size(); ++u){
                const uint64_t ku = cur[u].key;
                const long long du= cur[u].dist;
                int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCMaxT::digit_at(ku,i,powB,B);
                for(int i=0;i<m;++i){
                    int aii=ai[i];
                    if(aii<=0) continue;
                    if(i<m-1 && (aii-1) < ai[i+1]) continue; // 回退后仍非增
                    int jprev = aii - 1;

                    // 静态窗口存在性检查（落点必须允许是 Lcur）
                    if(LBwin && UBwin){
                        int lb = (*LBwin)[i][jprev];
                        int ub = (*UBwin)[i][jprev];
                        if(!(lb <= Lcur && Lcur <= ub)) continue;
                    }

                    uint64_t keyNew = lightOCMaxT::encode_digit_dec(ku,i,powB);
                    long long c = (long long)(Lcur) * weight_ij(i, jprev, m, h) * dV;
                    long long nd = du + c;

                    auto it = il.find(keyNew);
                    if(it==il.end()){
                        int v=(int)nl.size(); nl.push_back({keyNew, nd}); il.emplace(keyNew, v);
                    }else{
                        int v=it->second;
                        if(nd<nl[v].dist) nl[v].dist=nd;
                    }
                }
            }
        }
        for(size_t k=0;k<nxt_local.size();++k){
            for(auto &kvpair : nxt_local[k]){
                auto it = nxt_id.find(kvpair.key);
                if(it==nxt_id.end()){
                    int v=(int)next.size(); next.push_back(kvpair); nxt_id.emplace(kvpair.key, v);
                }else{
                    int v = it->second;
                    if(kvpair.dist<next[v].dist) next[v].dist=kvpair.dist;
                }
            }
            vector<KV>().swap(nxt_local[k]); unordered_map<uint64_t,int>().swap(id_local[k]);
        }
#else
        for(int u=0; u<(int)cur.size(); ++u){
            const uint64_t ku = cur[u].key;
            const long long du= cur[u].dist;
            int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCMaxT::digit_at(ku,i,powB,B);
            for(int i=0;i<m;++i){
                int aii=ai[i];
                if(aii<=0) continue;
                if(i<m-1 && (aii-1) < ai[i+1]) continue;
                int jprev = aii - 1;

                if(LBwin && UBwin){
                    int lb = (*LBwin)[i][jprev];
                    int ub = (*UBwin)[i][jprev];
                    if(!(lb <= Lcur && Lcur <= ub)) continue;
                }

                uint64_t keyNew = lightOCMaxT::encode_digit_dec(ku,i,powB);
                long long c = (long long)(Lcur) * weight_ij(i, jprev, m, h) * dV;
                long long nd = du + c;

                auto it = nxt_id.find(keyNew);
                if(it==nxt_id.end()){
                    int v=(int)next.size(); next.push_back({keyNew, nd}); nxt_id.emplace(keyNew, v);
                }else{
                    int v = it->second;
                    if(nd<next[v].dist) next[v].dist=nd;
                }
            }
        }
#endif
        if(progress && (vlevel>=2) && (Lcur%16==0)){
            LOG(2, [&](std::ofstream& o){
                o<<"[DP-bwd] level "<<Lcur<<" states="<<cur.size()<<" -> next="<<next.size()<<"\n";
            });
        }
        cur.swap(next);
    }
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv : cur) out.emplace(kv.key, kv.dist);
}

// ------- Hirschberg + mid-shift，带 cut 兼容性过滤 -------
static bool find_mid_intersection(
    int m,int h, long long dV, int Tcap,
    const vector<vector<int>>* LBwin, const vector<vector<int>>* UBwin,
    const vector<uint64_t>& powB, uint64_t B,
    bool progress, int omp_threads, int vlevel,
    int l_lo, uint64_t key_lo, int l_hi, uint64_t key_hi,
    int l_mid, uint64_t& out_key_mid)
{
    unordered_map<uint64_t,long long> F, Bmap;
    forward_to_level_map(m,h,key_lo,l_lo,l_mid-l_lo,dV,Tcap,LBwin,UBwin,powB,B,F,progress,omp_threads,vlevel);
    backward_to_level_map(m,h,key_hi,l_hi,l_hi-l_mid,dV,Tcap,LBwin,UBwin,powB,B,Bmap,progress,omp_threads,vlevel);

    long long best = std::numeric_limits<long long>::max();
    uint64_t key_mid = 0;

    auto try_key = [&](uint64_t key, long long d1, long long d2)->void{
        // cut 兼容性过滤
        vector<int> a(m);
        for(int i=0;i<m;++i) a[i] = lightOCMaxT::digit_at(key,i,powB,B);
        if(!cut_state_T_compatible(m,h,l_mid,a,Tcap,*LBwin,*UBwin)) return;
        long long v = d1 + d2;
        if(v<best){ best=v; key_mid=key; }
    };

    if(F.size()<=Bmap.size()){
        for(auto &kv : F){
            auto it = Bmap.find(kv.first);
            if(it==Bmap.end()) continue;
            try_key(kv.first, kv.second, it->second);
        }
    }else{
        for(auto &kv : Bmap){
            auto it = F.find(kv.first);
            if(it==F.end()) continue;
            try_key(kv.first, it->second, kv.second);
        }
    }
    if(best==std::numeric_limits<long long>::max()){
        LOG(2, [&](std::ofstream& o){
            o<<"[MID] empty at l_mid="<<l_mid<<" F="<<F.size()<<" B="<<Bmap.size()
             <<" (after cut-compat)\n";
        });
        return false;
    }
    out_key_mid = key_mid;
    LOG(3, [&](std::ofstream& o){
        o<<"[MID] hit at l_mid="<<l_mid<<" cost="<<best<<"  F="<<F.size()<<" B="<<Bmap.size()<<"\n";
    });
    return true;
}

static void reconstruct_hirschberg(
    int m,int h, long long dV, int Tcap,
    const vector<vector<int>>* LBwin,
    const vector<vector<int>>* UBwin,
    const vector<uint64_t>& powB, uint64_t B,
    bool progress, int omp_threads, int vlevel,
    int l_lo, uint64_t key_lo, int l_hi, uint64_t key_hi,
    vector<vector<int>>& y)
{
    int len = l_hi - l_lo;
    if(len==0) return;
    if(len==1){
        int sel_i=-1;
        for(int i=0;i<m;++i){
            int a_lo = lightOCMaxT::digit_at(key_lo,i,powB,B);
            int a_hi = lightOCMaxT::digit_at(key_hi,i,powB,B);
            if(a_hi == a_lo+1){ sel_i=i; break; }
        }
        if(sel_i<0) throw std::runtime_error("Base step mismatch");
        int j0 = lightOCMaxT::digit_at(key_lo, sel_i, powB, B);
        y[sel_i][j0] = l_lo + 1;
        return;
    }

    int l_mid = l_lo + len/2;
    uint64_t key_mid=0;
    if(!find_mid_intersection(m,h,dV,Tcap,LBwin,UBwin,powB,B,progress,omp_threads,vlevel,
                              l_lo,key_lo,l_hi,key_hi,l_mid,key_mid))
    {
        const int SHIFT_MAX = 8;
        bool ok=false;
        for(int d=1; d<=SHIFT_MAX && !ok; ++d){
            if(l_mid-d>l_lo){
                if(find_mid_intersection(m,h,dV,Tcap,LBwin,UBwin,powB,B,progress,omp_threads,vlevel,
                                         l_lo,key_lo,l_hi,key_hi,l_mid-d,key_mid)){
                    l_mid -= d; ok=true; break;
                }
            }
            if(l_mid+d<l_hi){
                if(find_mid_intersection(m,h,dV,Tcap,LBwin,UBwin,powB,B,progress,omp_threads,vlevel,
                                         l_lo,key_lo,l_hi,key_hi,l_mid+d,key_mid)){
                    l_mid += d; ok=true; break;
                }
            }
        }
        if(!ok){
            LOG(1, [&](std::ofstream& o){
                o<<"[MID] FAIL: no midpoint in ["<<l_lo<<","<<l_hi<<"], tried shift ±"<<SHIFT_MAX<<"\n";
            });
            throw std::runtime_error("No midpoint found");
        }
    }

    reconstruct_hirschberg(m,h,dV,Tcap,LBwin,UBwin,powB,B,progress,omp_threads,vlevel,
                           l_lo,key_lo,l_mid,key_mid,y);
    reconstruct_hirschberg(m,h,dV,Tcap,LBwin,UBwin,powB,B,progress,omp_threads,vlevel,
                           l_mid,key_mid,l_hi,key_hi,y);
}

// --------------- solve ------------------
OCMaxTResult lightOCMaxT::solve(int m, int h, int T_in){
    LOGOPEN(cfg_.logfile, cfg_.log_append, cfg_.verbose_level);

    const int n_full = m*h;
    int T = T_in;
    if(T<=0) T=n_full;
    if(T>n_full) T=n_full;

    const uint64_t B = (uint64_t)h + 1u;
    vector<uint64_t> powB(m);
    powB[m-1] = 1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t tmp = (__uint128_t)powB[k+1]*(__uint128_t)B;
        if(tmp > (__uint128_t)std::numeric_limits<uint64_t>::max())
            throw std::runtime_error("state key overflow");
        powB[k]=(uint64_t)tmp;
    }

    uint64_t key0=0ull, keyN=0ull;
    for(int i=0;i<m;++i) keyN += powB[i]*(uint64_t)h;

    LOG(1, [&](std::ofstream& o){
        o<<"[OC-SpanT] n="<<n_full<<" T="<<T
#ifdef _OPENMP
         <<" thr="<<(cfg_.omp_threads>0?cfg_.omp_threads:omp_get_max_threads())
#endif
         <<" v="<<cfg_.verbose_level<<" (Δ≤T)\n";
    });

    // 生成 Δ≤T 窗口
    vector<vector<int>> LBwin, UBwin;
    compute_windows_spanT(m,h,T,n_full,LBwin,UBwin,cfg_.verbose_level);

    // Hirschberg 分治（前向有动态 Δ 充分约束；后向无动态约束；mid 处做 cut 兼容性过滤）
    vector<vector<int>> y(m, vector<int>(h,0));
    reconstruct_hirschberg(m,h,cfg_.dV,T,&LBwin,&UBwin,powB,B,
                           cfg_.progress,cfg_.omp_threads,cfg_.verbose_level,
                           0,key0,n_full,keyN,y);

    // 评估与校验
    long long total = hpwl_full_neighbors(y, cfg_.dV);
    int dmax=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) dmax=std::max(dmax, y[i][j+1]-y[i][j]);
    for(int j=0;j<h;++j) for(int i=0;i+1<m;++i) dmax=std::max(dmax, y[i+1][j]-y[i][j]);
    if(dmax>T){
        LOG(1, [&](std::ofstream& o){ o<<"[POST] maxΔ="<<dmax<<" > T="<<T<<"\n"; });
        throw std::runtime_error("Post-check failed: max Δ > T");
    }
    LOG(1, [&](std::ofstream& o){ o<<"[POST] HPWL="<<total<<" maxΔ="<<dmax<<" OK\n"; });

    OCMaxTResult R;
    R.m=m; R.h=h; R.n_full=n_full; R.T=T; R.dV=cfg_.dV;
    R.y_order=std::move(y);
    R.hpwl_full=total;
    R.total_cost=total;
    R.hpwl_prefix=0;
    R.oc_ok=true;
    return R;
}
