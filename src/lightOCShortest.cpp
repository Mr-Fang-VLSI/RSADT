#include "lightOCShortest.h"
#include <unordered_map>
#include <iostream>
#include <limits>
#include <algorithm>
#include <stdexcept>

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
            throw std::runtime_error("State key overflows 64-bit; switch to vector<int> key for this m,h");
        }
        powB[k] = (uint64_t)tmp;
    }

    // 层容器（仅持久化各层用于回溯）
    vector<vector<Entry>> layers; layers.reserve(n+1);
    layers.emplace_back();           // level 0
    layers.back().push_back(Entry{0ull, 0ll, -1, -1});

    // 主循环：逐层 DP（仅读上一层，构造下一层）
    for(int lvl=0; lvl<n; ++lvl){
        const auto& cur = layers.back(); // 只读
        vector<Entry> next;              next.reserve(std::max<size_t>(cur.size()*2, 16));
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(std::max<size_t>(cur.size()*2, 16));

        for(int u=0; u<(int)cur.size(); ++u){
            const uint64_t key_u = cur[u].key;
            const long long du  = cur[u].dist;

            for(int i=0;i<m;++i){
                const int ai = digit_at(key_u, i, powB, B);
                if(ai >= h) continue;
                if(i>0){
                    const int aim1 = digit_at(key_u, i-1, powB, B);
                    if(ai + 1 > aim1) continue;
                }
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

        if (cfg_.progress && (lvl % 16 == 0)) {
            cout << "[DP] level " << lvl << " states=" << cur.size()
                 << " -> next=" << next.size() << endl;
        }

        if (next.empty()) {
            throw std::runtime_error("Empty next layer at level " + std::to_string(lvl));
        }
        // 持久化下一层（不要再 swap 回来！）
        layers.emplace_back(std::move(next));
    }

    // 选择最后一层的最短距离状态
    const auto& last = layers.back();
    int best_v = 0; long long best_d = last[0].dist;
    for(int v=1; v<(int)last.size(); ++v){
        if(last[v].dist < best_d){ best_d = last[v].dist; best_v = v; }
    }

    // 回溯构造 y_order（只回放 T 层；若 T<n_full，剩余为 0）
    vector<vector<int>> y(m, vector<int>(h, 0));
    int v = best_v;
    for(int lvl=n; lvl>=1; --lvl){
        const auto& L  = layers[lvl];
        const auto& P  = layers[lvl-1];
        const int u    = L[v].parent_idx;
        const int i    = L[v].move_i;
        const uint64_t key_u = P[u].key;
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
             << ", last_states=" << last.size()
             << ", total_cost=" << R.total_cost
             << ", HPWL=" << R.hpwl
             << ", OC=" << (R.oc_ok?"OK":"FAIL") << endl;
    }
    return R;
}
