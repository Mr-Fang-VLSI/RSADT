#include "ParamThresholdOCPlacer.h"
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <queue>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <limits>

using std::vector;
using std::cout;
using std::endl;

static inline int vid(int i, int j, int h) { return i*h + j; }

long long ParamThresholdOCPlacer::hpwl_edges_sum(const vector<vector<int>>& y,
                                                 const vector<long long>& x,
                                                 long long dV) {
    (void)x; // 单列 x=0
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long S = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j + 1 < h; ++j)
            S += (long long)std::llabs(y[i][j+1] - y[i][j]) * dV; // 水平邻边（dx=0）
    for (int i = 0; i + 1 < m; ++i)
        for (int j = 0; j < h; ++j)
            S += (long long)std::llabs(y[i+1][j] - y[i][j]) * dV; // 竖向邻边
    return S;
}

bool ParamThresholdOCPlacer::check_global_OC(const vector<vector<int>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
    for (int j1 = 0; j1 < h; ++j1)
    for (int i2 = i1; i2 < m; ++i2)
    for (int j2 = j1; j2 < h; ++j2)
        if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

// 在半整数阈值序列上做最大闭包（min-cut）
// lam2 ∈ {5,3,1,-1,-3} 对应 λ = {2.5,1.5,0.5,-0.5,-1.5}
// w2 = 2*w ∈ {-4,-2,0,2,4}
// s->v cap = max(0, w2 - lam2), v->t cap = max(0, lam2 - w2)
// 覆盖弧：v -> 前驱 INF
// 额外：强制 keep[v]=1 的点始终在源侧（加 s->v INF）
static void max_closure_lambda2(int m, int h,
                                const vector<int>& w2,
                                int lam2,
                                const vector<char>& keep, // 强制包含（保证嵌套）
                                vector<char>& inS,        // 输出：源侧闭包
                                long long& cuts_run,
                                bool verbose)
{
    const int n = m*h;
    const long long INF = 1e12;

    lemon::ListDigraph g;
    auto s = g.addNode();
    auto t = g.addNode();

    vector<lemon::ListDigraph::Node> V(n);
    for (int v = 0; v < n; ++v) V[v] = g.addNode();

    lemon::ListDigraph::ArcMap<long long> cap(g);
    auto add_arc = [&](lemon::ListDigraph::Node u, lemon::ListDigraph::Node v, long long c){
        auto a = g.addArc(u, v); cap[a] = c; return a;
    };

    long long edge_cnt = 0;
    long long sumPos = 0;

    for (int v = 0; v < n; ++v) {
        long long p = (long long)w2[v] - (long long)lam2; // 可能为负
        if (p >= 0) { add_arc(s, V[v], p); sumPos += p; ++edge_cnt; }
        else        { add_arc(V[v], t, -p); ++edge_cnt; }
        if (keep[v]) { add_arc(s, V[v], INF); ++edge_cnt; } // 强制留在源侧
    }

    // 覆盖弧（OC 前驱）
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int v = vid(i,j,h);
            if (i > 0) add_arc(V[v], V[vid(i-1,j,h)], INF), ++edge_cnt; // 向下闭合
            if (j > 0) add_arc(V[v], V[vid(i,j-1,h)], INF), ++edge_cnt;
        }

    if (verbose) {
        cout << "[λ-sweep] lam2=" << lam2
             << " nodes=" << (n+2)
             << " arcs~=" << edge_cnt << endl;
    }

    lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<long long>> pf(g, cap, s, t);
    pf.runMinCut(); ++cuts_run;

    inS.assign(n, 0);
    int sz = 0; long long sumVal = 0;
    for (int v = 0; v < n; ++v) {
        bool onS = pf.minCut(V[v]);
        inS[v] = onS ? 1 : 0;
        if (onS) {
            long long p = (long long)w2[v] - (long long)lam2;
            sumVal += p; ++sz;
        }
    }
    if (verbose) {
        cout << "  [cut] |S|=" << sz
             << " sumVal=" << sumVal
             << " sumPos(src)=" << sumPos << endl;
    }
}

PTResult ParamThresholdOCPlacer::solve(int m, int h) {
    const int n = m*h;
    // 2*权重
    vector<int> w2(n, 0);
    for (int i = 0, v=0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v)
            w2[v] = 2 * base_weight(i, j, m, h); // ∈{-4,-2,0,2,4}

    // 半整数阈值序列（降序）
    const int lam2_list[5] = {5, 3, 1, -1, -3};

    // 嵌套链保障：keep 集合单调扩张
    vector<char> keep(n, 0);
    vector<int> lstar2(n, std::numeric_limits<int>::min()); // 该点出现的最大 lam2
    long long cuts = 0;

    for (int k = 0; k < 5; ++k) {
        int lam2 = lam2_list[k];
        vector<char> inS; inS.reserve(n);
        max_closure_lambda2(m, h, w2, lam2, keep, inS, cuts, cfg_.verbose);

        // 更新 lstar2（首次出现时记录 lam2）
        for (int v = 0; v < n; ++v) {
            if (inS[v] && lstar2[v] == std::numeric_limits<int>::min())
                lstar2[v] = lam2;
        }
        // 维护嵌套：下一轮强制保留本轮 S
        keep = inS;
    }

    // 若某些点在所有阈值都未出现（极端内部 w=0 的点可能不在任何半整数上出现），
    // 给它们分配最小 lam2 作为后备，保证可比较性（不会影响最优性）。
    for (int v = 0; v < n; ++v)
        if (lstar2[v] == std::numeric_limits<int>::min())
            lstar2[v] = -3; // 最小阈值

    // —— 按 lstar2 降序做 Kahn 拓扑排序（严格满足 OC） —— //
    vector<int> indeg(n, 0);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int v = vid(i,j,h);
            if (i > 0) ++indeg[v];
            if (j > 0) ++indeg[v];
        }

    struct NodeKey { int lstar2; int i; int j; int v; };
    struct Cmp {
        bool operator()(const NodeKey& A, const NodeKey& B) const {
            if (A.lstar2 != B.lstar2) return A.lstar2 < B.lstar2; // 大的先出
            // 次关键字可偏向“更靠上/更靠右”的点先排（对 HPWL 更友好）
            if (A.i != B.i) return A.i < B.i; // i大者（靠上）优先 => 用 <
            return A.j < B.j;                 // j大者（靠右）优先 => 用 <
        }
    };
    std::priority_queue<NodeKey, std::vector<NodeKey>, Cmp> pq;

    for (int i = 0, v=0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v)
            if (indeg[v] == 0) pq.push(NodeKey{lstar2[v], i, j, v});

    vector<int> order; order.reserve(n);
    while (!pq.empty()) {
        NodeKey nk = pq.top(); pq.pop();
        order.push_back(nk.v);
        int i = nk.i, j = nk.j;
        // 后继：上(i+1,j)、右(i,j+1)
        if (i + 1 < m) {
            int u = vid(i+1, j, h);
            if (--indeg[u] == 0) pq.push(NodeKey{lstar2[u], i+1, j, u});
        }
        if (j + 1 < h) {
            int u = vid(i, j+1, h);
            if (--indeg[u] == 0) pq.push(NodeKey{lstar2[u], i, j+1, u});
        }
    }

    // 回填 y_order（按 order 顺序赋 rank）
    vector<vector<int>> Y(m, vector<int>(h, 0));
    for (int t = 0; t < (int)order.size(); ++t) {
        int v = order[t];
        int i = v / h, j = v % h;
        Y[i][j] = t + 1;
    }

    PTResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(Y);
    R.x_of_col.assign(h, 0);

    // 成本与核验
    long long Z = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int bw = base_weight(i,j,m,h);
            Z += (long long)bw * (long long)R.y_order[i][j] * cfg_.dV;
        }
    R.cost = Z;
    R.actual_hpwl = hpwl_edges_sum(R.y_order, R.x_of_col, cfg_.dV);
    R.oc_ok = check_global_OC(R.y_order);

    // 唯一占用 (x,y)
    {
        std::unordered_set<uint64_t> occ; occ.reserve((size_t)n*2);
        bool ok = true;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < h; ++j) {
                uint64_t x = 0ull;
                uint64_t yv = (uint64_t)R.y_order[i][j];
                uint64_t key64 = (x << 32) ^ yv;
                if (!occ.insert(key64).second) ok = false;
            }
        R.unique_ok = ok;
    }

    R.cuts_run = (int)cuts;

    if (cfg_.verbose) {
        cout << "[λ-sweep] cuts=" << R.cuts_run
             << ", cost=" << R.cost
             << ", HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok? "OK":"FAIL")
             << ", Unique=" << (R.unique_ok? "OK":"FAIL") << endl;
    }
    return R;
}
