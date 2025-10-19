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

struct Entry {
    uint64_t key;       // base-(h+1) 编码后的 a[0..m-1]
    long long dist;     // 到达该状态的最短距离
    int parent_idx;     // 上一层的 index
    int16_t move_i;     // 从 parent 到该状态时选择的行 i
};

OCShortestResult lightOCShortest::solve(int m, int h){
    const int n_full = m*h;
    const int n = (cfg_.max_steps>0 && cfg_.max_steps<n_full) ? cfg_.max_steps : n_full;

    // base 与幂表（powB[idx] = B^(m-1-idx)）
    const uint64_t B = (uint64_t)h + 1u;
    vector<uint64_t> powB(m);
    powB[m-1] = 1ull;
    for(int k=m-2;k>=0;--k){
        __uint128_t tmp = (__uint128_t)powB[k+1] * (__uint128_t)B;
        if (tmp > ( (__uint128_t)std::numeric_limits<uint64_t>::max() )) {
            throw std::runtime_error("State key overflows 64-bit; switch to vector<int> key");
        }
        powB[k] = (uint64_t)tmp;
    }

    // 层容器（仅持久化各层用于回溯）
    vector<vector<Entry>> layers; layers.reserve(n+1);
    layers.emplace_back();           // level 0
    layers.back().push_back(Entry{0ull, 0ll, -1, -1});

    // 线程配置
    int OMP_N = 1;
#ifdef _OPENMP
    OMP_N = (cfg_.omp_threads > 0 ? cfg_.omp_threads : omp_get_max_threads());
    if (cfg_.omp_threads > 0) omp_set_num_threads(cfg_.omp_threads);
#endif

    // 主循环：逐层 DP
    for(int lvl=0; lvl<n; ++lvl){
        const auto& cur = layers.back(); // 只读上一层
        if (cur.empty()) throw std::runtime_error("Empty current layer");

        if (cfg_.progress && (lvl % 16 == 0)) {
            cout << "[DP] level " << lvl << " states=" << cur.size();
#ifdef _OPENMP
            cout << " (threads="<< OMP_N <<")";
#endif
            cout << endl;
        }

        // --- 分片并行：每个线程本地去重 ---
        // 预估本层的扩展规模（粗略）：每个父状态大概 1 + #drops 个子状态
        size_t est = std::max<size_t>(cur.size()*4, 16);
#ifdef _OPENMP
        vector< vector<Entry> >                   nxt_local(OMP_N);
        vector< unordered_map<uint64_t,int> >     map_local(OMP_N);
        for (int t=0; t<OMP_N; ++t) {
            nxt_local[t].reserve(est / OMP_N + 16);
            map_local[t].reserve(est / OMP_N + 16);
        }

#pragma omp parallel
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            auto &next = nxt_local[tid];
            auto &nid  = map_local[tid];

#pragma omp for schedule(static)
            for (int u=0; u<(int)cur.size(); ++u){
                const uint64_t key_u = cur[u].key;
                const long long du  = cur[u].dist;

                // 一次性解码 a[0..m-1]
                int ai_buf[32]; // m<=32
                for (int i=0;i<m;++i) {
                    ai_buf[i] = digit_at(key_u, i, powB, B);
                }

                // 只在合法 i 上扩展：i==0 或 ai+1 <= a_{i-1}
                for (int i=0;i<m;++i){
                    const int ai = ai_buf[i];
                    if (ai >= h) continue;
                    if (i>0 && ai + 1 > ai_buf[i-1]) continue;

                    const uint64_t key_v = encode_digit_inc(key_u, i, powB);
                    const int j_new = ai;
                    const long long c = (long long)(lvl+1) * weight_ij(i, j_new, m, h) * cfg_.dV;
                    const long long nd = du + c;

                    auto it = nid.find(key_v);
                    if(it == nid.end()){
                        int v = (int)next.size();
                        next.push_back(Entry{ key_v, nd, u, (int16_t)i });
                        nid.emplace(key_v, v);
                    }else{
                        int v = it->second;
                        if(nd < next[v].dist){
                            next[v].dist = nd;
                            next[v].parent_idx = u;
                            next[v].move_i = (int16_t)i;
                        }
                    }
                }
            } // end for u
        } // end parallel

        // --- 归并各线程的 next → 全局 next（无锁阶段） ---
        vector<Entry> next; next.reserve(est);
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(est);

        for (int t=0; t<OMP_N; ++t) {
            auto &loc_map = map_local[t];
            auto &loc_vec = nxt_local[t];
            for (auto &kv : loc_map) {
                uint64_t key = kv.first;
                int li = kv.second;
                const Entry &E = loc_vec[li];

                auto it = nxt_id.find(key);
                if (it == nxt_id.end()) {
                    int gi = (int)next.size();
                    next.push_back(E);
                    nxt_id.emplace(key, gi);
                } else {
                    int gi = it->second;
                    if (E.dist < next[gi].dist) {
                        next[gi].dist = E.dist;
                        next[gi].parent_idx = E.parent_idx;
                        next[gi].move_i = E.move_i;
                        next[gi].key = E.key;
                    }
                }
            }
            // 释放线程本地容量（交给系统回收）
            vector<Entry>().swap(loc_vec);
            unordered_map<uint64_t,int>().swap(loc_map);
        }
#else
        // --- 单线程路径（无 OpenMP） ---
        vector<Entry> next; next.reserve(est);
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(est);

        for (int u=0; u<(int)cur.size(); ++u){
            const uint64_t key_u = cur[u].key;
            const long long du  = cur[u].dist;

            int ai_buf[32];
            for (int i=0;i<m;++i) ai_buf[i] = digit_at(key_u, i, powB, B);

            for (int i=0;i<m;++i){
                const int ai = ai_buf[i];
                if (ai >= h) continue;
                if (i>0 && ai + 1 > ai_buf[i-1]) continue;

                const uint64_t key_v = encode_digit_inc(key_u, i, powB);
                const int j_new = ai;
                const long long c = (long long)(lvl+1) * weight_ij(i, j_new, m, h) * cfg_.dV;
                const long long nd = du + c;

                auto it = nxt_id.find(key_v);
                if(it == nxt_id.end()){
                    int v = (int)next.size();
                    next.push_back(Entry{ key_v, nd, u, (int16_t)i });
                    nxt_id.emplace(key_v, v);
                }else{
                    int v = it->second;
                    if(nd < next[v].dist){
                        next[v].dist = nd;
                        next[v].parent_idx = u;
                        next[v].move_i = (int16_t)i;
                    }
                }
            }
        }
#endif

        if (next.empty()) {
            throw std::runtime_error("Empty next layer at level " + std::to_string(lvl));
        }
        layers.emplace_back(std::move(next));
    } // end for lvl

    // 选择最后一层的最短距离状态
    const auto& last = layers.back();
    int best_v = 0; long long best_d = last[0].dist;
    for(int v2=1; v2<(int)last.size(); ++v2){
        if(last[v2].dist < best_d){ best_d = last[v2].dist; best_v = v2; }
    }

    // 回溯构造 y_order
    vector<vector<int>> y(m, vector<int>(h, 0));
    int v = best_v;
    for(int lvl=n; lvl>=1; --lvl){
        const auto& L  = layers[lvl];
        const auto& P  = layers[lvl-1];
        const int u    = L[v].parent_idx;
        const int i    = L[v].move_i;
        const uint64_t key_u = P[u].key;
        // 扩展前列号 j0
        const int j0   = digit_at(key_u, i, powB, B);
        y[i][j0] = lvl;
        v = u;
    }

    OCShortestResult R;
    R.m=m; R.h=h; R.n=n_full; R.dV=cfg_.dV;
    R.y_order = std::move(y);
    R.total_cost = best_d;
    R.hpwl = hpwl_sum(R.y_order, cfg_.dV);
    R.oc_ok = check_OC(R.y_order);

    if(cfg_.verbose){
        cout << "[OC-Shortest] levels=" << n
#ifdef _OPENMP
             << ", threads=" << OMP_N
#endif
             << ", last_states=" << last.size()
             << ", total_cost=" << R.total_cost
             << ", HPWL=" << R.hpwl
             << ", OC=" << (R.oc_ok?"OK":"FAIL") << endl;
    }
    return R;
}
