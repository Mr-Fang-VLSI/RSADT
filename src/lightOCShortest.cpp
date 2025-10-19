#include "lightOCShortest.h"
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

long long lightOCShortest::hpwl_sum(const vector<vector<int>>& y, long long dV){
    const int m = (int)y.size(), h = (int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i)
        for(int j=0;j+1<h;++j) S += llabsll((long long)y[i][j+1]-y[i][j])*dV;
    for(int i=0;i+1<m;++i)
        for(int j=0;j<h;++j)   S += llabsll((long long)y[i+1][j]-y[i][j])*dV;
    return S;
}
bool lightOCShortest::check_OC(const vector<vector<int>>& y){
    const int m = (int)y.size(), h=(int)y[0].size();
    for(int i1=0;i1<m;++i1)
        for(int j1=0;j1<h;++j1)
            for(int i2=i1;i2<m;++i2)
                for(int j2=j1;j2<h;++j2)
                    if(y[i1][j1] > y[i2][j2]) return false;
    return true;
}

struct KV { uint64_t key; long long dist; };

// ===== 两层滚动：正向到目标层，返回该层 dist map =====
static void forward_to_level_map(
    int m,int h, uint64_t key0,int l0,int steps, long long dV,
    const vector<uint64_t>& powB, uint64_t B,
    unordered_map<uint64_t,long long>& out,
    bool progress, int omp_threads)
{
    vector<KV> cur; cur.reserve(4096);
    cur.push_back({key0, 0});
    for(int lvl=0; lvl<steps; ++lvl){
        // 并行展开到下一层
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
                for(int i=0;i<m;++i) ai[i] = lightOCShortest::digit_at(ku,i,powB,B);
                for(int i=0;i<m;++i){
                    int aii = ai[i];
                    if(aii>=h) continue;
                    if(i>0 && aii+1>ai[i-1]) continue;
                    uint64_t kv = lightOCShortest::encode_digit_inc(ku,i,powB);
                    long long c = (long long)(l0 + lvl + 1) * weight_ij(i, aii, m, h) * dV;
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
        // 合并线程本地
        for(size_t t=0;t<nxt_local.size();++t){
            for(auto &kv : nxt_local[t]){
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
            vector<KV>().swap(nxt_local[t]);
            unordered_map<uint64_t,int>().swap(id_local[t]);
        }
#else
        for(int u=0; u<(int)cur.size(); ++u){
            const uint64_t ku = cur[u].key;
            const long long du= cur[u].dist;
            int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCShortest::digit_at(ku,i,powB,B);
            for(int i=0;i<m;++i){
                int aii=ai[i];
                if(aii>=h) continue;
                if(i>0 && aii+1>ai[i-1]) continue;
                uint64_t kv = lightOCShortest::encode_digit_inc(ku,i,powB);
                long long c = (long long)(l0 + lvl + 1) * weight_ij(i, aii, m, h) * dV;
                long long nd = du + c;
                auto it = nxt_id.find(kv);
                if(it==nxt_id.end()){
                    int v=(int)next.size();
                    next.push_back({kv, nd});
                    nxt_id.emplace(kv, v);
                }else{
                    int v=it->second;
                    if(nd<next[v].dist) next[v].dist=nd;
                }
            }
        }
#endif
        if(progress && ((l0+lvl)%16==0))
            cout << "[DP-fwd] level " << (l0+lvl) << " states="<<cur.size()<<" -> next="<<next.size()<<endl;

        cur.swap(next);
    }
    // 输出哈希：key -> dist
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv : cur) out.emplace(kv.key, kv.dist);
}

// ===== 两层滚动：反向到目标层，返回该层 dist map =====
static void backward_to_level_map(
    int m,int h, uint64_t keyN,int lN,int steps, long long dV,
    const vector<uint64_t>& powB, uint64_t B,
    unordered_map<uint64_t,long long>& out,
    bool progress, int omp_threads)
{
    vector<KV> cur; cur.reserve(4096);
    cur.push_back({keyN, 0});
    for(int t=0; t<steps; ++t){
        int Lcur = lN - t; // 这一步的正向 rank = Lcur
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
                int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCShortest::digit_at(ku,i,powB,B);
                for(int i=0;i<m;++i){
                    int aii=ai[i];
                    if(aii<=0) continue;
                    if(i<m-1 && (aii-1) < ai[i+1]) continue; // 保持非增
                    uint64_t kv = lightOCShortest::encode_digit_dec(ku,i,powB);
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
            int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCShortest::digit_at(ku,i,powB,B);
            for(int i=0;i<m;++i){
                int aii=ai[i];
                if(aii<=0) continue;
                if(i<m-1 && (aii-1) < ai[i+1]) continue;
                uint64_t kv = lightOCShortest::encode_digit_dec(ku,i,powB);
                int jprev = aii - 1;
                long long c = (long long)(Lcur) * weight_ij(i, jprev, m, h) * dV;
                long long nd = du + c;
                auto it = nxt_id.find(kv);
                if(it==nxt_id.end()){
                    int v=(int)next.size();
                    next.push_back({kv, nd});
                    nxt_id.emplace(kv, v);
                }else{
                    int v=it->second;
                    if(nd<next[v].dist) next[v].dist=nd;
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

// ===== 递归回溯：Hirschberg =====
static void reconstruct_hirschberg(
    int m,int h, long long dV,
    const vector<uint64_t>& powB, uint64_t B, bool progress, int omp_threads,
    int l_lo, uint64_t key_lo, int l_hi, uint64_t key_hi,
    vector<vector<int>>& y)
{
    int len = l_hi - l_lo;
    if(len==0) return;
    if(len==1){
        // 找出唯一变化的行
        int sel_i = -1;
        for(int i=0;i<m;++i){
            int a_lo = lightOCShortest::digit_at(key_lo,i,powB,B);
            int a_hi = lightOCShortest::digit_at(key_hi,i,powB,B);
            if(a_hi == a_lo + 1){ sel_i = i; break; }
        }
        if(sel_i<0) throw std::runtime_error("Base step mismatch");
        int j0 = lightOCShortest::digit_at(key_lo, sel_i, powB, B);
        y[sel_i][j0] = l_lo + 1;
        return;
    }

    int l_mid = l_lo + len/2;
    // F: 0..mid ； B: hi..mid
    unordered_map<uint64_t,long long> F, Bmap;
    forward_to_level_map(m,h, key_lo, l_lo, l_mid-l_lo, dV, powB, B, F, progress, omp_threads);
    backward_to_level_map(m,h, key_hi, l_hi, l_hi-l_mid, dV, powB, B, Bmap, progress, omp_threads);

    // 选择中点 key*
    uint64_t key_mid = 0; long long best = std::numeric_limits<long long>::max();
    // 让小的表驱动遍历
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

    // 递归左右半段
    reconstruct_hirschberg(m,h,dV,powB,B,progress,omp_threads, l_lo, key_lo, l_mid, key_mid, y);
    reconstruct_hirschberg(m,h,dV,powB,B,progress,omp_threads, l_mid, key_mid, l_hi, key_hi, y);
}

OCShortestResult lightOCShortest::solve(int m, int h){
    const int n_full = m*h;
    const int n = (cfg_.max_steps>0 && cfg_.max_steps<n_full) ? cfg_.max_steps : n_full;

    // base & pow table
    const uint64_t B = (uint64_t)h + 1u;
    vector<uint64_t> powB(m);
    powB[m-1] = 1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t tmp = (__uint128_t)powB[k+1] * (__uint128_t)B;
        if(tmp > (__uint128_t)std::numeric_limits<uint64_t>::max())
            throw std::runtime_error("State key overflows 64-bit; use vector<int> encoding");
        powB[k] = (uint64_t)tmp;
    }

    // 起点/终点编码
    uint64_t key0 = 0ull;
    uint64_t keyN = 0ull; // 所有 a[i]=h
    for(int i=0;i<m;++i) keyN += powB[i] * (uint64_t)h;

    vector<vector<int>> y(m, vector<int>(h, 0));

    if(cfg_.low_mem){
        // Hirschberg 低内存回溯
        if(cfg_.verbose){
            cout << "[OC-Shortest/low-mem] m="<<m<<" h="<<h<<" levels="<<n
#ifdef _OPENMP
                 << " threads=" << (cfg_.omp_threads>0?cfg_.omp_threads:omp_get_max_threads())
#endif
                 << endl;
        }
        reconstruct_hirschberg(m,h,cfg_.dV,powB,B,cfg_.progress,cfg_.omp_threads, 0,key0, n,keyN, y);
        // 求总成本：用端点化简（或邻接求和），两者在 OC 下等价
        long long total = hpwl_sum(y, cfg_.dV);
        OCShortestResult R{m,h,n,cfg_.dV, std::move(y), total, total, check_OC(y)};
        if(cfg_.verbose){
            cout << "[Verify] total_cost="<<R.total_cost<<", HPWL="<<R.hpwl
                 << ", OC="<<(R.oc_ok?"OK":"FAIL")<<endl;
        }
        return R;
    }
    else{
        // （保持你原“全层持久化”版，略。建议 16×16 时启用 low_mem）
        throw std::runtime_error("non-low_mem path not implemented in this build");
    }
}
