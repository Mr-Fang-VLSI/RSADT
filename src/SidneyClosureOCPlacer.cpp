#include "SidneyClosureOCPlacer.h"
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <limits>

using std::vector;
using std::cout;
using std::endl;

static inline int node_id(int i, int j, int h) { return i*h + j; }

// 真值 HPWL
long long SidneyClosureOCPlacer::hpwl_edges_sum(const vector<vector<int>>& y,
                                                const vector<long long>& x,
                                                long long dV) {
    (void)x; // 单列
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long S = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j+1 < h; ++j)
            S += (long long)std::llabs(y[i][j+1]-y[i][j]) * dV;
    for (int i = 0; i+1 < m; ++i)
        for (int j = 0; j < h; ++j)
            S += (long long)std::llabs(y[i+1][j]-y[i][j]) * dV;
    return S;
}

bool SidneyClosureOCPlacer::check_global_OC(const vector<vector<int>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
    for (int j1 = 0; j1 < h; ++j1)
    for (int i2 = i1; i2 < m; ++i2)
    for (int j2 = j1; j2 < h; ++j2)
        if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

// —— 用一次最小割求“给定 λ 的最大闭包” ——
// 返回：inS（闭包 S），sumVal = ∑_{v∈S}(W[v] - λ)
static std::pair<int, long long>
max_closure_at_lambda(int m, int h,
                      const vector<char>& alive,
                      const vector<long long>& W,
                      long long lambda,
                      long long INF,
                      int &flows_counter,
                      vector<char> &inS_out)
{
    const int n = m*h;

    lemon::ListDigraph g;
    auto s = g.addNode();
    auto t = g.addNode();
    vector<lemon::ListDigraph::Node> V(n, lemon::INVALID);

    for (int v = 0; v < n; ++v) if (alive[v]) V[v] = g.addNode();

    lemon::ListDigraph::ArcMap<long long> cap(g);
    auto add_arc = [&](lemon::ListDigraph::Node u, lemon::ListDigraph::Node v, long long c){
        auto a = g.addArc(u, v); cap[a] = c; return a;
    };

    long long sumPos = 0;
    for (int i = 0, v=0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v) if (alive[v]) {
            long long w = W[v] - lambda;
            if (w >= 0) { add_arc(s, V[v], w); sumPos += w; }
            else        { add_arc(V[v], t, -w); }
        }

    // 前驱约束：若选 (i,j) 则必须选其前驱（下、左）
    for (int i = 0, v=0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v) if (alive[v]) {
            if (i > 0) {
                int pre = node_id(i-1, j, h);
                if (alive[pre]) add_arc(V[v], V[pre], INF);
            }
            if (j > 0) {
                int pre = node_id(i, j-1, h);
                if (alive[pre]) add_arc(V[v], V[pre], INF);
            }
        }

    // TODO: maxT 约束：按需要在此加入额外 INF 边（单列）

    lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<long long>> pf(g, cap, s, t);
    pf.runMinCut(); ++flows_counter;

    vector<char> inS(n, 0);
    int sz = 0; long long sumVal = 0;
    for (int v = 0; v < n; ++v) if (alive[v]) {
        bool onS = pf.minCut(V[v]);
        if (onS) { inS[v] = 1; ++sz; sumVal += (W[v] - lambda); }
    }
    inS_out.swap(inS);
    return {sz, sumVal};
}

// —— 二分 λ 找最大密度闭包（非空），返回 inS ——
// 令 BASE 远大于 n，用 W[v] = base_weight*BASE + eps[v] 保证词典序稳定。
// 令 predicate(mid): “存在非空闭包 S 使 ∑(W - mid) ≥ 0”。
// 取“最大使 predicate 为真”的 mid，并返回其闭包（非空）。
static vector<char>
max_density_closure(int m, int h,
                    const vector<char>& alive,
                    const vector<long long>& W,
                    long long &flows_counter)
{
    const int n = m*h;
    const long long INF = 1e12;

    long long lo = std::numeric_limits<long long>::max();
    long long hi = std::numeric_limits<long long>::min();
    for (int v = 0; v < n; ++v) if (alive[v]) {
        lo = std::min(lo, W[v]);
        hi = std::max(hi, W[v]);
    }
    // 放宽边界，确保：predicate(lo-1) 为真（至少全体闭包），predicate(hi+1) 为假（空集）
    lo = lo - 1;
    hi = hi + 1;

    vector<char> bestS; bestS.reserve(n);
    // 二分：取上中位，找最大的 mid
    while (lo < hi) {
        long long mid = lo + (hi - lo + 1)/2;
        vector<char> inS; auto [sz, val] = max_closure_at_lambda(m,h,alive,W,mid,INF,(int&)flows_counter,inS);

        bool ok = (sz > 0) && (val >= 0); // 非空且值≥0
        if (ok) { lo = mid; bestS = std::move(inS); }
        else    { hi = mid - 1; }
    }

    if (bestS.empty()) {
        // 极端场景：全负且任何密度都<lo；取 lo-1 再求一次得到非空闭包（最“不差”的）
        long long mid = lo - 1;
        vector<char> inS; auto [sz, val] = max_closure_at_lambda(m,h,alive,W,mid,INF,(int&)flows_counter,inS);
        if (sz == 0) {
            // 仍然空（理论上不应发生），强制取任意可行最小点（比如左下角）
            vector<char> fallback(n,0);
            for (int v=0; v<n; ++v) if (alive[v]) { fallback[v]=1; break; }
            return fallback;
        }
        return inS;
    }
    return bestS;
}

SidneyResult SidneyClosureOCPlacer::solve(int m, int h) {
    const int n = m*h;
    const long long BASE = n + 7; // 词典序放大（主目标优先）

    // 构造带扰动的 W
    vector<long long> W(n,0);
    for (int i = 0, v=0; i<m; ++i)
        for (int j = 0; j<h; ++j, ++v) {
            int bw = base_weight(i,j,m,h);
            long long eps = (long long)v; // 唯一化
            W[v] = (long long)bw*BASE + eps;
        }

    // 递归 Sidney 分解：每次取最大密度闭包放前面
    vector<char> alive(n, 1);
    vector<int> order; order.reserve(n);
    long long flows = 0;

    // 用 lambda 封装递归
    std::function<void(const vector<char>&)> dfs =
    [&](const vector<char>& A){
        // 统计剩余
        int rest = 0; for (char c : A) if (c) ++rest;
        if (rest == 0) return;
        // 取最大密度闭包（非空）
        vector<char> S = max_density_closure(m,h,A,W,flows);

        // 拆分 S 与补
        vector<char> B_in = S;            // 子图 S
        vector<char> B_out = A;           // 子图补
        for (int v=0; v<n; ++v) if (A[v]) {
            if (S[v]) B_out[v]=0; else B_in[v]=0;
        }
        // 先递归 S，再递归补
        dfs(B_in);
        dfs(B_out);
    };

    // 把 dfs 改成“记录秩”的形式：在叶子处返回单点
    // 为了显式构造线性序，我们写个包装器：
    std::function<void(const vector<char>&)> build =
    [&](const vector<char>& A){
        int rest = 0; for (char c : A) if (c) ++rest;
        if (rest == 0) return;
        if (rest == 1) {
            for (int v=0; v<n; ++v) if (A[v]) { order.push_back(v); break; }
            return;
        }
        vector<char> S = max_density_closure(m,h,A,W,flows);
        vector<char> A_in = S, A_out = A;
        for (int v=0; v<n; ++v) if (A[v]) {
            if (S[v]) A_out[v]=0; else A_in[v]=0;
        }
        build(A_in);
        build(A_out);
    };

    build(alive);

    // 回填 y_order（按 order 顺序赋 rank）
    vector<vector<int>> Y(m, vector<int>(h,0));
    for (int t=0; t<(int)order.size(); ++t) {
        int v = order[t];
        int i = v / h, j = v % h;
        Y[i][j] = t + 1;
    }

    SidneyResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(Y);
    R.x_of_col.assign(h, 0);

    // 统计 cost 与核验
    long long Z = 0;
    for (int i=0, v=0; i<m; ++i)
        for (int j=0; j<h; ++j, ++v){
            int bw = base_weight(i,j,m,h);
            Z += (long long)bw * (long long) ( /*rank=*/R.y_order[i][j] ) * cfg_.dV;
        }
    R.cost = Z;
    R.actual_hpwl = hpwl_edges_sum(R.y_order, R.x_of_col, cfg_.dV);
    R.oc_ok = check_global_OC(R.y_order);

    // 唯一占用
    {
        std::unordered_set<uint64_t> occ; occ.reserve((size_t)n*2);
        bool ok = true;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < h; ++j) {
                uint64_t x=0, yv=(uint64_t)R.y_order[i][j];
                uint64_t key = (x<<32)^yv;
                if (!occ.insert(key).second) ok=false;
            }
        R.unique_ok = ok;
    }

    R.flows_run = (int)flows;

    if (cfg_.verbose) {
        cout << "[Sidney-Closure] flows=" << R.flows_run
             << ", cost=" << R.cost
             << ", HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok? "OK":"FAIL")
             << ", Unique=" << (R.unique_ok? "OK":"FAIL") << endl;
    }
    return R;
}
