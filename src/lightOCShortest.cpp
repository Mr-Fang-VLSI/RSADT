#include "lightOCShortest.h"
#include <unordered_map>
#include <iostream>
#include <limits>
#include <algorithm>

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

// 每层持久化的条目：仅存压缩状态key、到达距离、父在上一层的位置、扩展行
struct Entry {
    uint64_t key;          // base-(h+1) 编码后的 a[0..m-1]
    long long dist;        // 到达该状态的最短距离
    int parent_idx;        // 上一层的 index
    int16_t move_i;        // 从 parent 到该状态时选择的行 i
};

OCShortestResult lightOCShortest::solve(int m, int h){
    const int n_full = m*h;
    const int n = (cfg_.max_steps>0 && cfg_.max_steps<n_full) ? cfg_.max_steps : n_full;

    // base 与幂表
    const uint64_t B = (uint64_t)h + 1u;
    vector<uint64_t> powB(m);
    // powB[idx] = B^(m-1-idx)
    powB[m-1] = 1ull;
    for(int k=m-2;k>=0;--k) powB[k] = powB[k+1] * B;

    // 层容器（持久化，用于回溯），只存本层条目
    vector<vector<Entry>> layers; layers.reserve(n+1);

    // 当前层 & 下一层（去重用 map）
    vector<Entry> cur; cur.reserve(4096);
    vector<Entry> nxt; nxt.reserve(16384);
    unordered_map<uint64_t,int> cur_id, nxt_id;
    cur_id.reserve(4096); nxt_id.reserve(16384);

    // 初始层（全 0，key=0）
    cur.push_back(Entry{0ull, 0ll, -1, -1});
    cur_id.emplace(0ull, 0);
    layers.push_back(cur);

    // 主循环：逐层 DP
    for(int lvl=0; lvl<n; ++lvl){
        nxt.clear(); nxt_id.clear();

        for(int u=0; u<(int)cur.size(); ++u){
            const uint64_t key_u = cur[u].key;
            const long long du  = cur[u].dist;

            // 预取 digit（一次除法展开全部行）
            // 为减少除法次数，这里只在需要的位置取 digit
            // 检查每一行是否可扩（a[i]<h 且 a[i]+1<=a[i-1]）
            for(int i=0;i<m;++i){
                const int ai   = digit_at(key_u, i, powB, B);
                if(ai >= h) continue;
                if(i>0){
                    const int aim1 = digit_at(key_u, i-1, powB, B);
                    if(ai + 1 > aim1) continue;
                }
                // 生成子状态
                const uint64_t key_v = encode_digit_inc(key_u, i, powB);
                const int j_new = ai; // 扩展前列号
                const long long c = (long long)(lvl+1) * weight_ij(i, j_new, m, h) * cfg_.dV;
                const long long nd = du + c;

                auto it = nxt_id.find(key_v);
                if(it == nxt_id.end()){
                    int v = (int)nxt.size();
                    nxt.push_back(Entry{ key_v, nd, u, (int16_t)i });
                    nxt_id.emplace(key_v, v);
                }else{
                    int v = it->second;
                    if(nd < nxt[v].dist){
                        nxt[v].dist = nd;
                        nxt[v].parent_idx = u;
                        nxt[v].move_i = (int16_t)i;
                    }
                }
            }
        }

        // 进度
        if(cfg_.progress && (lvl%16==0)){
            cout << "[DP] level " << lvl << " states=" << cur.size() << " -> next=" << nxt.size() << endl;
        }

        layers.push_back(std::move(nxt));            // 保存当前“下一层”
        nxt = vector<Entry>(); nxt.reserve(16384);   // 复位容器，避免深拷贝
        cur.swap(layers.back());                     // cur 指向刚保存的层
        cur_id.swap(nxt_id);                         // map 也对调
    }

    // 从最后一层选最小 dist 的状态
    const auto& last = layers.back();
    if(last.empty()){
        throw std::runtime_error("Empty last layer (no feasible state)");
    }
    int best_v = 0; long long best_d = last[0].dist;
    for(int v=1; v<(int)last.size(); ++v){
        if(last[v].dist < best_d){ best_d = last[v].dist; best_v = v; }
    }

    // 回溯构造 y_order
    vector<vector<int>> y(m, vector<int>(h, 0));
    int v = best_v;
    for(int lvl=n; lvl>=1; --lvl){
        const auto& L  = layers[lvl];
        const auto& Lp = layers[lvl-1];
        const int u = L[v].parent_idx;
        const int i = L[v].move_i;
        const uint64_t key_u = Lp[u].key;
        const int j0 = digit_at(key_u, i, powB, B); // 扩展前列号
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
