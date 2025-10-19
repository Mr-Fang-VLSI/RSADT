#include "AssignOCMcmfPlacer.h"
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <cassert>

using std::vector;
using std::cout;
using std::endl;

long long AssignOCMcmfPlacer::hpwl_edges_sum(const vector<vector<int>>& y,
                                             long long dV) {
    int m = (int)y.size(), h = (int)y[0].size();
    long long S = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j + 1 < h; ++j)
            S += (long long)std::llabs(y[i][j+1] - y[i][j]) * dV;
    for (int i = 0; i + 1 < m; ++i)
        for (int j = 0; j < h; ++j)
            S += (long long)std::llabs(y[i+1][j] - y[i][j]) * dV;
    return S;
}

bool AssignOCMcmfPlacer::check_global_OC(const vector<vector<int>>& y) {
    int m = (int)y.size(), h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
        for (int j1 = 0; j1 < h; ++j1)
            for (int i2 = i1; i2 < m; ++i2)
                for (int j2 = j1; j2 < h; ++j2)
                    if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

// 构建“前缀闭包网络”：层 r=1..n，每层 n 个点节点 V_{v,r}。
// R_r(供给1)->V_{v,r}(cost=w(v)*dV), 层内闭包 V_{v,r}->V_{u,r} (INF) 对覆盖边 u<-v，
// 层间单调 V_{v,r}->V_{v,r+1} (INF)，最顶层 V_{v,n}->t (cap=1)。
AssignOCResult AssignOCMcmfPlacer::solve(int m, int h) {
    using Dig = lemon::ListDigraph;
    Dig g;

    const int n = m*h;
    const long long INF = 1e9; // 足够大即可（边数小，费用小）
    auto s = g.addNode();
    auto t = g.addNode();

    // R_r 节点
    vector<Dig::Node> Rnode(n+1);
    for (int r = 1; r <= n; ++r) Rnode[r] = g.addNode();



    // V_{v,r} 节点
    vector<vector<Dig::Node>> V(n, vector<Dig::Node>(n+1));
    for (int v = 0; v < n; ++v)
        for (int r = 1; r <= n; ++r)
            V[v][r] = g.addNode();

    lemon::ListDigraph::ArcMap<int> cap(g);
    lemon::ListDigraph::ArcMap<long long> cost(g);
    lemon::ListDigraph::NodeMap<long long> supply(g);

    for (Dig::NodeIt it(g); it != lemon::INVALID; ++it) supply[it] = 0;

    auto add_arc = [&](Dig::Node u, Dig::Node v, int c, long long w)->Dig::Arc{
        auto a = g.addArc(u, v);
        cap[a] = c;
        cost[a] = w;
        return a;
    };

    long long arcs_cnt = 0;

    // 源→每层供给 R_r：1 单位
    for (int r = 1; r <= n; ++r) {
        add_arc(s, Rnode[r], 1, 0); ++arcs_cnt;
    }
    // 最顶层 V_{v,n}→汇：cap=1
    for (int v = 0; v < n; ++v) {
        add_arc(V[v][n], t, 1, 0); ++arcs_cnt;
    }

    // R_r→V_{v,r}（选择 v 加入 S_r）：cap=1, cost=w(v)*dV
    for (int i = 0, v = 0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v) {
            long long w = (long long)base_weight(i,j,m,h) * (long long)cfg_.dV;
            for (int r = 1; r <= n; ++r) {
                add_arc(Rnode[r], V[v][r], 1, w);
                ++arcs_cnt;
            }
        }

    // 层内闭包（覆盖边：上/左）
    auto add_cover = [&](int u, int v, int r){ // u<-v（v依赖u）
        add_arc(V[v][r], V[u][r], (int)INF, 0);
        ++arcs_cnt;
    };
    for (int i = 0, v = 0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v) {
            int u_up = (i>0 ? vid(i-1,j,h) : -1);
            int u_left = (j>0 ? vid(i,j-1,h) : -1);
            for (int r = 1; r <= n; ++r) {
                if (u_up   >= 0) add_cover(u_up,   v, r);
                if (u_left >= 0) add_cover(u_left, v, r);
            }
        }

    // 层间单调：V_{v,r}→V_{v,r+1} (INF)
    for (int v = 0; v < n; ++v)
        for (int r = 1; r < n; ++r) {
            add_arc(V[v][r], V[v][r+1], (int)INF, 0);
            ++arcs_cnt;
        }

    // supply：源 +n，汇 -n
    supply[s] =  n;
    supply[t] = -n;

    if (cfg_.verbose) {
        cout << "[build] nodes=" << lemon::countNodes(g)
             << ", arcs~=" << arcs_cnt
             << " (n=" << n << ", m=" << m << ", h=" << h << ")\n";
    }

    // 运行 NetworkSimplex
    lemon::NetworkSimplex<Dig, int, long long> ns(g);
    ns.upperMap(cap).costMap(cost).supplyMap(supply);
    auto st = ns.run();
    if (st != decltype(ns)::OPTIMAL) {
        throw std::runtime_error("NetworkSimplex did not find OPTIMAL.");
    }
    long long min_cost = ns.totalCost();

    // 复原：每层 r 选了哪个 v ∈ S_r？（看 R_r→V_{v,r} 的流）
    // 然后按 “第一次进入的层 r == y(v)” 回填秩。
    vector<int> first_layer(n, 0);
    for (int r = 1; r <= n; ++r) {
        for (lemon::ListDigraph::OutArcIt a(g, Rnode[r]); a != lemon::INVALID; ++a) {
            int flow = ns.flow(a);
            if (flow == 1) {
                // a: R[r] -> V[v][r]
                Dig::Node v_node = g.target(a);
                // 反查 v
                int v_id = -1;
                // 用 i/j 二重循环反查（开销可接受；n<=~1000 时影响极小）
                for (int v = 0; v < n; ++v) {
                    if (V[v][r] == v_node) { v_id = v; break; }
                }
                if (v_id < 0) continue;
                first_layer[v_id] = (first_layer[v_id] == 0 ? r : std::min(first_layer[v_id], r));
            }
        }
    }

    vector<vector<int>> Y(m, vector<int>(h, 0));
    // first_layer[v] 就是 y(v)
    for (int v = 0; v < n; ++v) {
        int i = v / h, j = v % h;
        int y = first_layer[v];
        if (y <= 0) {
            // 极端不该发生：说明该点从未被任何层选择（与供给相悖）
            y = n; 
        }
        Y[i][j] = y;
    }

    AssignOCResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(Y);
    R.x_of_col.assign(h, 0);

    long long Z = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int w = base_weight(i,j,m,h);
            Z += (long long)w * (long long)R.y_order[i][j] * cfg_.dV;
        }
    R.cost = Z;
    R.actual_hpwl = hpwl_edges_sum(R.y_order, cfg_.dV);
    R.oc_ok = check_global_OC(R.y_order);

    // 唯一占用（y 唯一）
    {
        std::unordered_set<int> occ; occ.reserve((size_t)n*2);
        bool ok = true;
        for (int i=0;i<m;++i) for (int j=0;j<h;++j) {
            int yv = R.y_order[i][j];
            if (!occ.insert(yv).second) ok=false;
        }
        R.unique_ok = ok;
    }

    R.nodes = lemon::countNodes(g);
    R.arcs  = arcs_cnt;

    if (cfg_.verbose) {
        cout << "[assign+OC MCMF] nodes=" << R.nodes
             << ", arcs~=" << R.arcs
             << ", cost=" << R.cost
             << ", HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok? "OK":"FAIL")
             << ", Unique=" << (R.unique_ok? "OK":"FAIL") << endl;
    }
    return R;
}
