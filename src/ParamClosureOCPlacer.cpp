#include "ParamClosureOCPlacer.h"
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

long long ParamClosureOCPlacer::hpwl_edges_sum(const vector<vector<int>>& y,
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

bool ParamClosureOCPlacer::check_global_OC(const vector<vector<int>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
    for (int j1 = 0; j1 < h; ++j1)
    for (int i2 = i1; i2 < m; ++i2)
    for (int j2 = j1; j2 < h; ++j2)
        if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

PCResult ParamClosureOCPlacer::solve(int m, int h) {
    const int n = m * h;
    const long long BASE = n + 7;      // 词典序放大系数
    const long long INF  = 1e12;       // min-cut 的“无限容量”
    // λ 搜索范围（覆盖 base_weight*BASE 与 eps 的总幅度）
    const long long LAM_MAX =  (2LL*BASE + n + 10);
    const long long LAM_MIN = -(2LL*BASE + n + 10);

    // 构图（一次性）：顶点=格点；边=“若选 v 则必须选其前驱”
    lemon::ListDigraph g;
    lemon::ListDigraph::Node s = g.addNode();
    lemon::ListDigraph::Node t = g.addNode();

    vector<lemon::ListDigraph::Node> V(n);
    for (int i = 0; i < n; ++i) V[i] = g.addNode();

    // 记录 Source/Sink 边，以便每次 λ 变更时只改 capacity
    vector<lemon::ListDigraph::Arc> A_src(n), A_snk(n);

    lemon::ListDigraph::ArcMap<long long> cap(g);
    auto add_arc = [&](lemon::ListDigraph::Node u, lemon::ListDigraph::Node v, long long c){
        auto a = g.addArc(u, v);
        cap[a] = c;
        return a;
    };

    // 节点 id 与权重（带极小扰动的词典序放大）
    vector<long long> W(n, 0);
    {
        int idx = 0;
        for (int ii = 0; ii < m; ++ii)
            for (int jj = 0; jj < h; ++jj, ++idx) {
                int bw = base_weight(ii, jj, m, h);     // -2..+2
                long long eps = (long long)idx;         // 唯一化小扰动
                W[idx] = (long long)bw * BASE + eps;
                // 先放一条“空容量”的源/汇边（后续按 λ 改 cap）
                A_src[idx] = add_arc(s, V[idx], 0);
                A_snk[idx] = add_arc(V[idx], t, 0);
            }
    }

    // 前驱约束边：若选 (i,j) 则必须选 (i-1,j) 与 (i,j-1)
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int u = node_id(i,j,h);
            if (i > 0) {
                int pre = node_id(i-1, j, h);
                add_arc(V[u], V[pre], INF); // 选 u ⇒ 必选 pre
            }
            if (j > 0) {
                int pre = node_id(i, j-1, h);
                add_arc(V[u], V[pre], INF);
            }
        }

    // 闭包求解器（给定 λ，返回 S(λ) 及其大小）
    auto closure_at_lambda = [&](long long lam, vector<char>& inS) -> int {
        // 更新源/汇边容量：正接源，负接汇
        for (int id = 0; id < n; ++id) {
            long long w = W[id] - lam;
            if (w >= 0) {
                cap[A_src[id]] = w;
                cap[A_snk[id]] = 0;
            } else {
                cap[A_src[id]] = 0;
                cap[A_snk[id]] = -w;
            }
        }
        // max-flow / min-cut
        lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<long long>> pf(g, cap, s, t);
        pf.runMinCut(); // 对 Preflow，runMinCut() 就够；或 pf.run() 也可
        inS.assign(n, 0);
        int sz = 0;
        for (int id = 0; id < n; ++id) {
            bool onS = pf.minCut(V[id]); // true 表示在 source 一侧
            inS[id] = onS ? 1 : 0;
            if (onS) ++sz;
        }
        return sz;
    };

    // 递归分解：输出严格嵌套链，得到每个点的 rank（1..n）
    vector<int> rank(n, 0);
    int next_rank = 1;
    int flows = 0;

    // 为了保证二分能推进，我们先算两端闭包
    vector<char> S_hi, S_lo;
    int s_hi = closure_at_lambda(LAM_MAX, S_hi); ++flows; // 预期 0
    int s_lo = closure_at_lambda(LAM_MIN, S_lo); ++flows; // 预期 n

    auto assign_between = [&](auto&& self,
                              long long lam_hi, const vector<char>& S_hi, int s_hi,
                              long long lam_lo, const vector<char>& S_lo, int s_lo) -> void {
        if (s_lo == s_hi) return;
        if (s_lo == s_hi + 1) {
            // 差集里唯一的那个点拿到 rank = next_rank
            int id = -1;
            for (int v = 0; v < n; ++v) if (S_lo[v] && !S_hi[v]) { id = v; break; }
            rank[id] = next_rank++;
            return;
        }
        // 取中 λ；若闭包大小没有进展，就微调 λ（最多试几次）
        long long lam_mid = (lam_hi + lam_lo) / 2;
        vector<char> S_mid; int s_mid = -1;
        const int MAX_ADJ = 4*n + 10; // 足够小的尝试范围
        int tries = 0;
        while (tries < MAX_ADJ) {
            s_mid = closure_at_lambda(lam_mid, S_mid); ++flows;
            if (s_mid > s_hi && s_mid < s_lo) break;
            // 若等于边界，向内收缩
            if (s_mid <= s_hi) { lam_mid = lam_mid - 1; }
            else               { lam_mid = lam_mid + 1; }
            ++tries;
        }
        if (s_mid <= s_hi || s_mid >= s_lo) {
            // 仍然没有切入中间（极少发生），直接线性找一个切入值
            long long lam = lam_lo - 1;
            for (;;) {
                s_mid = closure_at_lambda(lam, S_mid); ++flows;
                if (s_mid > s_hi && s_mid < s_lo) { lam_mid = lam; break; }
                lam -= 1;
                if (lam < LAM_MIN - 10*BASE) {
                    throw std::runtime_error("[ParamClosure] cannot find mid lambda.");
                }
            }
        }
        // 递归左右两段
        self(self, lam_hi, S_hi, s_hi, lam_mid, S_mid, s_mid);
        self(self, lam_mid, S_mid, s_mid, lam_lo, S_lo, s_lo);
    };

    assign_between(assign_between, LAM_MAX, S_hi, s_hi, LAM_MIN, S_lo, s_lo);

    // 回填 y_order
    vector<vector<int>> Y(m, vector<int>(h, 0));
    for (int i = 0, id = 0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++id)
            Y[i][j] = rank[id];

    PCResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(Y);
    R.x_of_col.assign(h, 0);

    // 成本与核验
    long long Z = 0;
    for (int i = 0, id = 0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++id) {
            int bw = base_weight(i,j,m,h);
            Z += (long long)bw * (long long)rank[id] * cfg_.dV;
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

    R.flows_run = flows;

    if (cfg_.verbose) {
        cout << "[ParamClosure] flows=" << R.flows_run
             << ", cost=" << R.cost
             << ", HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok? "OK":"FAIL")
             << ", Unique=" << (R.unique_ok? "OK":"FAIL") << endl;
    }
    return R;
}
