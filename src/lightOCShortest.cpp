#include "lightOCShortest.h"
#include <unordered_map>
#include <queue>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <algorithm>

using std::vector;
using std::pair;
using std::cout;
using std::endl;

static inline long long llabsll(long long x){ return x>=0?x:-x; }

// 目标的权重（与你之前一致：四边-1/+1）
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
    long long S = 0;
    for(int i=0;i<m;++i)
        for(int j=0;j+1<h;++j)
            S += llabsll((long long)y[i][j+1] - y[i][j]) * dV;
    for(int i=0;i+1<m;++i)
        for(int j=0;j<h;++j)
            S += llabsll((long long)y[i+1][j] - y[i][j]) * dV;
    return S;
}

bool lightOCShortest::check_OC(const vector<vector<int>>& y){
    const int m = (int)y.size(), h = (int)y[0].size();
    for(int i1=0;i1<m;++i1)
        for(int j1=0;j1<h;++j1)
            for(int i2=i1;i2<m;++i2)
                for(int j2=j1;j2<h;++j2)
                    if(y[i1][j1] > y[i2][j2]) return false;
    return true;
}

// 自定义 hash：vector<int> 作为 key，避免 uint64 溢出
struct VecHash {
    size_t operator()(const vector<int>& a) const noexcept {
        // 64-bit 混合（splitmix 风格）
        uint64_t h = 1469598103934665603ull;
        for(int v : a){
            uint64_t x = (uint64_t)(uint32_t)v + 0x9e3779b97f4a7c15ull;
            h ^= x + (h<<6) + (h>>2);
        }
        return (size_t)h;
    }
};

OCShortestResult lightOCShortest::solve(int m, int h){
    const int n = m * h;

    // 起始状态 a[i]=0（自底向上），非增：a[0]>=a[1]>=...>=a[m-1]
    vector<int> a0(m, 0);

    // 状态池
    struct Node {
        vector<int> a; // 非增阶梯
        int level;     // sum(a)
        int parent = -1;
        int move_i  = -1; // 从parent到该状态时选择的行 i
    };
    vector<Node> nodes; nodes.reserve(1'500'000);

    // 索引：state -> id
    std::unordered_map<vector<int>, int, VecHash> id_of;
    id_of.reserve(1'500'000);

    auto add_state = [&](const vector<int>& a, int lvl, int parent, int move_i)->int{
        auto it = id_of.find(a);
        if(it != id_of.end()){
            int id = it->second;
            // 若该 id 尚未设置 parent（首次访问）可以保留更短 parent，但 Dijkstra 会处理距离
            return id;
        }
        int id = (int)nodes.size();
        nodes.push_back(Node{a, lvl, parent, move_i});
        id_of.emplace(nodes.back().a, id);
        return id;
    };

    auto sum_vec = [&](const vector<int>& a)->int{
        long long s=0; for(int v:a) s+=v; return (int)s;
    };

    int id_start = add_state(a0, 0, -1, -1);

    // Dijkstra
    const long long INF = std::numeric_limits<long long>::max()/4;
    vector<long long> dist; dist.reserve(1'500'000);
    vector<char>      seen; seen.reserve(1'500'000);

    // NOTE: nodes.size() 会动态增长，dist/seen 要与之同步扩张
    dist.resize(nodes.size(), INF);
    seen.resize(nodes.size(), 0);
    dist[id_start] = 0;

    struct QItem { long long d; int id; };
    auto cmp = [](const QItem& a, const QItem& b){ return a.d > b.d; };
    std::priority_queue<QItem, vector<QItem>, decltype(cmp)> pq(cmp);
    pq.push({0, id_start});

    int goal_id = -1;
    vector<int> b; b.reserve(m);

    // 为了避免重复计算，准备每个 i 的“是否可扩”快速检查
    auto legal_expand = [&](const vector<int>& a, int i)->bool{
        if(a[i] >= h) return false;
        if(i>0 && a[i] + 1 > a[i-1]) return false;
        return true;
    };

    // 主循环：搜索至 level==n 的任意状态
    while(!pq.empty()){
        auto [cd, u] = pq.top(); pq.pop();
        if(u >= (int)nodes.size()) continue; // 防御
        if(seen[u]) continue;
        if(cd != dist[u]) continue;
        seen[u] = 1;

        // 终点：level==n
        if(nodes[u].level == n){ goal_id = u; break; }

        const vector<int>& a = nodes[u].a;
        const int lvl = nodes[u].level; // 下一步秩 = lvl+1

        // 枚举可扩展行 i
        for(int i=0;i<m;++i){
            if(!legal_expand(a, i)) continue;
            b = a;
            int j_new = b[i]; // 扩展前列号
            b[i] = j_new + 1;
            // b 仍需保持非增，但由于只 +1，且前置检查 a[i]+1<=a[i-1] 已保证
            // 计算增量代价
            long long w = weight_ij(i, j_new, m, h);
            long long c = (long long)(lvl + 1) * w * cfg_.dV;

            int v = add_state(b, lvl+1, u, i);

            // 确保 dist/seen 可达
            if(v >= (int)dist.size()){
                dist.resize(v+1, INF);
                seen.resize(v+1, 0);
            }

            long long nd = cd + c;
            if(nd < dist[v]){
                dist[v] = nd;
                // 更新 parent 信息（保存在节点结构里）
                nodes[v].parent = u;
                nodes[v].move_i = i;
                pq.push({nd, v});
            }
        }
    }

    if(goal_id < 0){
        throw std::runtime_error("No goal state found (this should not happen).");
    }

    // 复原路径 → y_order
    vector<vector<int>> y(m, vector<int>(h, 0));
    // 从终点回溯到起点，沿途知道每一步选择了哪一行 i，对应 j = a[parent][i]（扩展前）
    vector<int> path_ids;
    for(int cur = goal_id; cur != -1; cur = nodes[cur].parent) path_ids.push_back(cur);
    std::reverse(path_ids.begin(), path_ids.end()); // 从起点到终点

    // 逐步回放：第 r 步（1-based）在 nodes[path_ids[r]] 相对 nodes[path_ids[r-1]] 的差异行 i
    for(int step = 1; step < (int)path_ids.size(); ++step){
        int v = path_ids[step];
        int u = path_ids[step-1];
        int i = nodes[v].move_i;
        int j0 = nodes[u].a[i];   // 扩展前的列号
        int r  = step;            // rank = step（1..n）
        y[i][j0] = r;
    }

    OCShortestResult R;
    R.m=m; R.h=h; R.n=n; R.dV=cfg_.dV;
    R.y_order = std::move(y);
    R.total_cost = dist[goal_id];
    R.hpwl = hpwl_sum(R.y_order, cfg_.dV);
    R.oc_ok = check_OC(R.y_order);

    if(cfg_.verbose){
        cout << "[OC-Shortest] nodes_explored=" << nodes.size()
             << ", total_cost=" << R.total_cost
             << ", HPWL=" << R.hpwl
             << ", OC=" << (R.oc_ok?"OK":"FAIL") << endl;
    }
    return R;
}
