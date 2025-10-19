#include "KClosureParamOCPlacer.h"
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <unordered_set>
#include <queue>
#include <iostream>
#include <algorithm>
#include <limits>
#include <cassert>

using std::vector;
using std::cout;
using std::endl;

long long KClosureParamOCPlacer::hpwl_edges_sum(const vector<vector<int>>& y,
                                                long long dV) {
    int m = (int)y.size(), h = (int)y[0].size();
    long long S = 0;
    // 水平邻边
    for (int i = 0; i < m; ++i)
        for (int j = 0; j + 1 < h; ++j)
            S += (long long)std::llabs(y[i][j+1] - y[i][j]) * dV;
    // 垂直邻边
    for (int i = 0; i + 1 < m; ++i)
        for (int j = 0; j < h; ++j)
            S += (long long)std::llabs(y[i+1][j] - y[i][j]) * dV;
    return S;
}

bool KClosureParamOCPlacer::check_global_OC(const vector<vector<int>>& y) {
    int m = (int)y.size(), h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
        for (int j1 = 0; j1 < h; ++j1)
            for (int i2 = i1; i2 < m; ++i2)
                for (int j2 = j1; j2 < h; ++j2)
                    if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

// 一次闭包（min-cut）：容量均为整数，keep=1 的点强制在源侧
int KClosureParamOCPlacer::max_closure_with_keep(
    int m, int h,
    const vector<long long>& W,
    long long lam,
    const vector<char>& keep,
    vector<char>& inS,
    int& cuts_counter,
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

    long long edge_cnt = 0, sumPos = 0;
    for (int v = 0; v < n; ++v) {
        long long p = W[v] - lam; // 可能为负
        if (p >= 0) { add_arc(s, V[v], p); sumPos += p; ++edge_cnt; }
        else        { add_arc(V[v], t, -p); ++edge_cnt; }
        if (keep[v]) { add_arc(s, V[v], INF); ++edge_cnt; }
    }
    // 覆盖弧（OC：选 v 必选其下、左前驱）
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int v = id(i,j,h);
            if (i > 0) add_arc(V[v], V[id(i-1,j,h)], INF), ++edge_cnt;
            if (j > 0) add_arc(V[v], V[id(i,j-1,h)], INF), ++edge_cnt;
        }

    if (verbose) {
        cout << "[closure] lam=" << lam
             << " nodes=" << (n+2)
             << " arcs~=" << edge_cnt
             << " keep=" << std::count(keep.begin(), keep.end(), (char)1) << endl;
    }

    lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<long long>> pf(g, cap, s, t);
    pf.runMinCut(); ++cuts_counter;

    inS.assign(n, 0);
    int sz = 0;
    for (int v = 0; v < n; ++v) {
        bool onS = pf.minCut(V[v]);
        inS[v] = onS ? 1 : 0;
        if (onS) ++sz;
    }
    if (verbose) cout << "  [cut] |S|=" << sz << " sumPos=" << sumPos << endl;
    return sz;
}

KcParamResult KClosureParamOCPlacer::solve(int m, int h) {
    const int n = m*h;

    // W = w * dV（整型）
    vector<long long> W(n, 0);
    for (int i = 0, v=0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v)
            W[v] = (long long)base_weight(i,j,m,h) * cfg_.dV;

    // λ 搜索边界：覆盖所有可能断点
    long long Wmax = 2 * cfg_.dV, Wmin = -2 * cfg_.dV;
    long long lam_lo = Wmin - 10*cfg_.dV; // 保守余量
    long long lam_hi = Wmax + 10*cfg_.dV;

    vector<char> keep(n, 0);     // 已选（S_t）
    vector<int>   rank(n, 0);    // 最终秩
    int next_rank = 1;
    int cuts = 0;

    // 逐 t=1..n
    for (int t = 1; t <= n; ++t) {
        // 二分 λ：找最大 λ 使 |S(λ)| >= t
        long long lo = lam_lo, hi = lam_hi;
        vector<char> bestS(n, 0); int bestSize = -1;

        for (int it = 0; it < cfg_.max_bisect_iter && lo < hi; ++it) {
            long long mid = lo + (hi - lo + 1)/2;
            vector<char> inS;
            int sz = max_closure_with_keep(m,h,W,mid,keep,inS,cuts,cfg_.verbose);
            if (sz >= t) { lo = mid; bestS = std::move(inS); bestSize = sz; }
            else         { hi = mid - 1; }
        }
        // 保险：取 lo 作为最终 λ
        vector<char> Slo; int szlo = max_closure_with_keep(m,h,W,lo,keep,Slo,cuts,cfg_.verbose);
        if (szlo < t) {
            // 极端边界（理论不应发生）：向下放宽
            vector<char> Sforce; (void)max_closure_with_keep(m,h,W,lo-1,keep,Sforce,cuts,cfg_.verbose);
            Slo.swap(Sforce); szlo = (int)std::count(Slo.begin(),Slo.end(),(char)1);
        }

        // 若 |Slo| == t，直接确定新增点
        // 若 |Slo| >  t，说明在该 λ 存在平台：在 Slo 内按“可加前沿”逐个补到 t
        vector<char>& St = Slo;
        int cur_sz = (int)std::count(keep.begin(), keep.end(), (char)1);
        assert(cur_sz == t-1);

        auto preds_ready = [&](int v)->bool{
            int i = v / h, j = v % h;
            if (i>0 && !keep[id(i-1,j,h)]) return false;
            if (j>0 && !keep[id(i,j-1,h)]) return false;
            return true;
        };

        // 从 St\keep 中挑选“所有前驱已在 keep”的点，逐个填到 t
        while (cur_sz < t) {
            int pick = -1;
            // 简单 tie-break：偏向“越靠上、越靠右”的点（经验上更利于 HPWL）
            int best_i = -1, best_j = -1;
            for (int v = 0; v < n; ++v) if (St[v] && !keep[v] && preds_ready(v)) {
                int i = v / h, j = v % h;
                if (i > best_i || (i == best_i && j > best_j)) {
                    best_i = i; best_j = j; pick = v;
                }
            }
            if (pick < 0) {
                // 不应发生：若发生，说明 St 不是 down-set，或者 keep 不可扩。回退：放宽 λ 再取一次。
                vector<char> Sfix; int szfix = max_closure_with_keep(m,h,W,lo-1,keep,Sfix,cuts,cfg_.verbose);
                (void)szfix; St = Sfix; continue;
            }
            keep[pick] = 1; rank[pick] = next_rank++; ++cur_sz;
        }
    }

    // 回填 y_order
    vector<vector<int>> Y(m, vector<int>(h,0));
    for (int v = 0; v < n; ++v) {
        int i = v / h, j = v % h;
        Y[i][j] = rank[v];
    }

    KcParamResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(Y);
    R.x_of_col.assign(h, 0);

    // 目标与校验
    long long Z = 0;
    for (int i=0;i<m;++i)
        for (int j=0;j<h;++j){
            int bw = base_weight(i,j,m,h);
            Z += (long long)bw * (long long)R.y_order[i][j] * cfg_.dV;
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
    R.cuts_total = cuts;

    if (cfg_.verbose) {
        cout << "[k-closure param] cuts=" << R.cuts_total
             << ", cost="    << R.cost
             << ", HPWL="    << R.actual_hpwl
             << ", OC="      << (R.oc_ok? "OK":"FAIL")
             << ", Unique="  << (R.unique_ok? "OK":"FAIL") << endl;
    }
    return R;
}
