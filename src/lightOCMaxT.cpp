#include "lightOCMaxT.h"
#include <unordered_map>
#include <iostream>
#include <limits>
#include <algorithm>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;
using std::unordered_map;
using std::cout;
using std::endl;

static inline long long weight_ij(int i,int j,int m,int h){
    long long w=0;
    if(i==0)     w -= 1;
    if(i==m-1)   w += 1;
    if(j==0)     w -= 1;
    if(j==h-1)   w += 1;
    return w;
}

static inline long long llabsll(long long x){ return x>=0?x:-x; }

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
bool lightOCMaxT::check_OC(const vector<vector<int>>& y){
    const int m = (int)y.size(), h=(int)y[0].size();
    for(int i1=0;i1<m;++i1)
        for(int j1=0;j1<h;++j1)
            for(int i2=i1;i2<m;++i2)
                for(int j2=j1;j2<h;++j2)
                    if(y[i1][j1] > y[i2][j2]) return false;
    return true;
}

struct KV { uint64_t key; long long dist; };

// 正向：从 (l0,key0) 走 steps 步，返回目标层的 key->dist
static void forward_to_level_map(
    int m,int h, uint64_t key0,int l0,int steps, long long dV,
    const vector<uint64_t>& powB, uint64_t B,
    unordered_map<uint64_t,long long>& out,
    bool progress, int omp_threads)
{
    vector<KV> cur; cur.reserve(4096);
    cur.push_back({key0, 0});
    for(int lvl=0; lvl<steps; ++lvl){
        size_t est = std::max<size_t>(cur.size()*4, 16);
        vector<KV> next; next.reserve(est);
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(est);

#ifdef _OPENMP
        int TH = omp_threads>0? omp_threads : omp_get_max_threads();
        vector< vector<KV> >                 nxt_local(TH);
        vector< unordered_map<uint64_t,int> >id_local(TH);
        for(int t=0;t<TH;++t){
            nxt_local[t].reserve(est/TH + 16);
            id_local[t].reserve(est/TH + 16);
        }
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
                int ai[32];
                for(int i=0;i<m;++i) ai[i] = lightOCMaxT::digit_at(ku,i,powB,B);
                for(int i=0;i<m;++i){
                    int aii = ai[i];
                    if(aii>=h) continue;
                    if(i>0 && aii+1>ai[i-1]) continue;
                    uint64_t kv = lightOCMaxT::encode_digit_inc(ku,i,powB);
                    long long c = (long long)(l0 + lvl + 1) * weight_ij(i, aii, m, h) * dV;
                    long long nd = du + c;
                    auto it = il.find(kv);
                    if(it==il.end()){
                        int v=(int)nl.size();
                        nl.push_back({kv, nd});
                        il.emplace(kv, v);
                    }else{
                        int v = it->second;
                        if(nd<nl[v].dist) nl[v].dist=nd;
                    }
                }
            }
        }
        for(size_t t=0;t<nxt_local.size();++t){
            for(auto &kv : nxt_local[t]){
                auto it = nxt_id.find(kv.key);
                if(it==nxt_id.end()){
                    int v=(int)next.size();
                    next.push_back(kv);
                    nxt_id.emplace(kv.key, v);
                }else{
                    int v = it->second;
                    if(kv.dist<next[v].dist) next[v].dist=kv.dist;
                }
            }
            vector<KV>().swap(nxt_local[t]);
            unordered_map<uint64_t,int>().swap(id_local[t]);
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
                uint64_t kv = lightOCMaxT::encode_digit_inc(ku,i,powB);
                long long c = (long long)(l0 + lvl + 1) * weight_ij(i, aii, m, h) * dV;
                long long nd = du + c;
                auto it = nxt_id.find(kv);
                if(it==nxt_id.end()){
                    int v=(int)next.size();
                    next.push_back({kv, nd});
                    nxt_id.emplace(kv, v);
                }else{
                    int v = it->second;
                    if(nd<next[v].dist) next[v].dist=nd;
                }
            }
        }
#endif
        if(progress && ((l0+lvl)%16==0))
            cout << "[DP-fwd] level " << (l0+lvl) << " states="<<cur.size()<<" -> next="<<next.size()<<endl;

        cur.swap(next);
    }
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv : cur) out.emplace(kv.key, kv.dist);
}

// 反向：从 (lN,keyN) 倒走 steps 步（合法前驱），返回目标层的 key->dist
static void backward_to_level_map(
    int m,int h, uint64_t keyN,int lN,int steps, long long dV,
    const vector<uint64_t>& powB, uint64_t B,
    unordered_map<uint64_t,long long>& out,
    bool progress, int omp_threads)
{
    vector<KV> cur; cur.reserve(4096);
    cur.push_back({keyN, 0});
    for(int t=0; t<steps; ++t){
        int Lcur = lN - t;
        size_t est = std::max<size_t>(cur.size()*4, 16);
        vector<KV> next; next.reserve(est);
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(est);

#ifdef _OPENMP
        int TH = omp_threads>0? omp_threads : omp_get_max_threads();
        vector< vector<KV> >                 nxt_local(TH);
        vector< unordered_map<uint64_t,int> >id_local(TH);
        for(int k=0;k<TH;++k){
            nxt_local[k].reserve(est/TH + 16);
            id_local[k].reserve(est/TH + 16);
        }
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
                    if(i<m-1 && (aii-1) < ai[i+1]) continue; // 保持非增
                    uint64_t kv = lightOCMaxT::encode_digit_dec(ku,i,powB);
                    int jprev = aii - 1;
                    long long c = (long long)(Lcur) * weight_ij(i, jprev, m, h) * dV;
                    long long nd = du + c;
                    auto it = il.find(kv);
                    if(it==il.end()){
                        int v=(int)nl.size();
                        nl.push_back({kv, nd});
                        il.emplace(kv, v);
                    }else{
                        int v=it->second;
                        if(nd<nl[v].dist) nl[v].dist=nd;
                    }
                }
            }
        }
        for(size_t k=0;k<nxt_local.size();++k){
            for(auto &kv : nxt_local[k]){
                auto it = nxt_id.find(kv.key);
                if(it==nxt_id.end()){
                    int v=(int)next.size();
                    next.push_back(kv);
                    nxt_id.emplace(kv.key, v);
                }else{
                    int v=it->second;
                    if(kv.dist<next[v].dist) next[v].dist=kv.dist;
                }
            }
            vector<KV>().swap(nxt_local[k]);
            unordered_map<uint64_t,int>().swap(id_local[k]);
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
                uint64_t kv = lightOCMaxT::encode_digit_dec(ku,i,powB);
                int jprev = aii - 1;
                long long c = (long long)(Lcur) * weight_ij(i, jprev, m, h) * dV;
                long long nd = du + c;
                auto it = nxt_id.find(kv);
                if(it==nxt_id.end()){
                    int v=(int)next.size();
                    next.push_back({kv, nd});
                    nxt_id.emplace(kv.key, v);
                }else{
                    int v=it->second;
                    if(nd<next[v].dist) next[v].dist=kv.dist;
                }
            }
        }
#endif
        if(progress && ((lN-t)%16==0))
            cout << "[DP-bwd] level " << (lN-t) << " states="<<cur.size()<<" -> next="<<next.size()<<endl;

        cur.swap(next);
    }
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv : cur) out.emplace(kv.key, kv.dist);
}

// 递归 Hirschberg 回溯：已知 (l_lo,key_lo) 与 (l_hi,key_hi)，写回 y 的 [l_lo+1..l_hi] 段
static void reconstruct_hirschberg(
    int m,int h, long long dV,
    const vector<uint64_t>& powB, uint64_t B, bool progress, int omp_threads,
    int l_lo, uint64_t key_lo, int l_hi, uint64_t key_hi,
    vector<vector<int>>& y)
{
    int len = l_hi - l_lo;
    if(len==0) return;
    if(len==1){
        int sel_i = -1;
        for(int i=0;i<m;++i){
            int a_lo = lightOCMaxT::digit_at(key_lo,i,powB,B);
            int a_hi = lightOCMaxT::digit_at(key_hi,i,powB,B);
            if(a_hi == a_lo + 1){ sel_i = i; break; }
        }
        if(sel_i<0) throw std::runtime_error("Base step mismatch");
        int j0 = lightOCMaxT::digit_at(key_lo, sel_i, powB, B);
        y[sel_i][j0] = l_lo + 1;
        return;
    }

    int l_mid = l_lo + len/2;
    unordered_map<uint64_t,long long> F, Bmap;
    forward_to_level_map(m,h, key_lo, l_lo, l_mid-l_lo, dV, powB, B, F, progress, omp_threads);
    backward_to_level_map(m,h, key_hi, l_hi, l_hi-l_mid, dV, powB, B, Bmap, progress, omp_threads);

    uint64_t key_mid = 0; long long best = std::numeric_limits<long long>::max();
    if(F.size() <= Bmap.size()){
        for(auto &kv : F){
            auto it = Bmap.find(kv.first);
            if(it==Bmap.end()) continue;
            long long v = kv.second + it->second;
            if(v < best){ best = v; key_mid = kv.first; }
        }
    }else{
        for(auto &kv : Bmap){
            auto it = F.find(kv.first);
            if(it==F.end()) continue;
            long long v = kv.second + it->second;
            if(v < best){ best = v; key_mid = kv.first; }
        }
    }
    if(best == std::numeric_limits<long long>::max())
        throw std::runtime_error("No midpoint found on shortest path");

    reconstruct_hirschberg(m,h,dV,powB,B,progress,omp_threads, l_lo, key_lo, l_mid, key_mid, y);
    reconstruct_hirschberg(m,h,dV,powB,B,progress,omp_threads, l_mid, key_mid, l_hi, key_hi, y);
}

OCMaxTResult lightOCMaxT::solve(int m, int h, int T_in){
    const int n_full = m*h;
    int T = T_in;
    if(T <= 0) T = n_full;
    if(T > n_full) T = n_full;

    const uint64_t B = (uint64_t)h + 1u;
    vector<uint64_t> powB(m);
    powB[m-1] = 1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t tmp = (__uint128_t)powB[k+1] * (__uint128_t)B;
        if(tmp > (__uint128_t)std::numeric_limits<uint64_t>::max())
            throw std::runtime_error("State key overflows 64-bit; use vector<int> encoding");
        powB[k] = (uint64_t)tmp;
    }

    uint64_t key0 = 0ull;
    uint64_t keyN = 0ull; for(int i=0;i<m;++i) keyN += powB[i] * (uint64_t)h;

    vector<vector<int>> y(m, vector<int>(h, 0));

    // 1) 如果 T < n_full：先找 T 层的“最佳终点” key_T^*
    uint64_t keyT = keyN; long long bestT = 0;
    if(T < n_full){
        unordered_map<uint64_t,long long> FT;
        forward_to_level_map(m,h, key0, 0, T, cfg_.dV, powB, B, FT, cfg_.progress, cfg_.omp_threads);

        keyT = 0; bestT = std::numeric_limits<long long>::max();
        for(auto &kv: FT){
            if(kv.second < bestT){ bestT = kv.second; keyT = kv.first; }
        }
        if(keyT==0 && bestT==std::numeric_limits<long long>::max())
            throw std::runtime_error("No feasible state at level T");
        // 2) 用 Hirschberg 在 [0..T] 上回溯
        reconstruct_hirschberg(m,h,cfg_.dV,powB,B,cfg_.progress,cfg_.omp_threads, 0,key0, T,keyT, y);
        // 3) 汇总
        OCMaxTResult R;
        R.m=m; R.h=h; R.n_full=n_full; R.T=T; R.dV=cfg_.dV;
        R.y_order = std::move(y);
        R.total_cost = bestT;
        R.hpwl_prefix = hpwl_prefix_neighbors(R.y_order, cfg_.dV);
        R.hpwl_full = 0; // 无意义（未完整赋值）
        R.oc_ok = check_OC(R.y_order);
        if(cfg_.verbose){
            cout << "[OC-MaxT] m="<<m<<" h="<<h<<" T="<<T
#ifdef _OPENMP
                 << " thr="<<(cfg_.omp_threads>0?cfg_.omp_threads:omp_get_max_threads())
#endif
                 << " total_cost="<<R.total_cost
                 << " HPWL_prefix="<<R.hpwl_prefix
                 << " OC="<<(R.oc_ok?"OK":"FAIL") << endl;
        }
        return R;
    }

    // 4) 否则（T==n_full）：回到完整求解（0..n_full）
    //    直接做一次 [0..n_full] 的 Hirschberg 回溯
    reconstruct_hirschberg(m,h,cfg_.dV,powB,B,cfg_.progress,cfg_.omp_threads, 0,key0, n_full,keyN, y);
    long long total = hpwl_full_neighbors(y, cfg_.dV);

    OCMaxTResult R;
    R.m=m; R.h=h; R.n_full=n_full; R.T=T; R.dV=cfg_.dV;
    R.y_order = std::move(y);
    R.total_cost = total;
    R.hpwl_full = total;
    R.hpwl_prefix = 0;
    R.oc_ok = check_OC(R.y_order);

    if(cfg_.verbose){
        cout << "[OC-Shortest] levels=" << n_full
#ifdef _OPENMP
             << " thr="<<(cfg_.omp_threads>0?cfg_.omp_threads:omp_get_max_threads())
#endif
             << " total_cost="<<R.total_cost
             << " HPWL="<<R.hpwl_full
             << " OC="<<(R.oc_ok?"OK":"FAIL") << endl;
    }
    return R;
}
