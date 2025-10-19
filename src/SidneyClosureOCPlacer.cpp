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
// 新增 verbose 打印与严格的 alive/V[v] 检查
static std::pair<int, long long>
max_closure_at_lambda(int m, int h,
                      const std::vector<char>& alive,
                      const std::vector<long long>& W,
                      long long lambda,
                      long long INF,
                      int &flows_counter,
                      std::vector<char> &inS_out,
                      bool verbose)
{
    const int n = m * h;
    // 统计 alive 顶点数量
    int alive_cnt = 0;
    for (int v = 0; v < n; ++v)
        if (alive[v]) ++alive_cnt;

    if (verbose) {
        std::cout << "[max_closure] lambda=" << lambda
                  << " |alive|=" << alive_cnt << std::endl;
    }

    if (alive_cnt <= 0) {
        inS_out.assign(n, 0);
        return {0, 0};
    }

    // 安全索引函数，避免越界
    auto safe_id = [&](int i, int j)->int {
        if (i < 0 || j < 0 || i >= m || j >= h)
            return -1;
        return i * h + j;
    };

    // -------- 以下保持原逻辑（构图/求 min-cut） --------
    lemon::ListDigraph g;
    lemon::ListDigraph::Node s = g.addNode();
    lemon::ListDigraph::Node t = g.addNode();

    std::vector<lemon::ListDigraph::Node> V(n, lemon::INVALID);
    for (int v = 0; v < n; ++v)
        if (alive[v]) V[v] = g.addNode();

    lemon::ListDigraph::ArcMap<long long> cap(g);
    auto add_arc = [&](lemon::ListDigraph::Node u,
                       lemon::ListDigraph::Node v,
                       long long c) {
        auto a = g.addArc(u, v);
        cap[a] = c;
        return a;
    };

    long long sumPos = 0;
    long long edge_count = 0;

    for (int v = 0; v < n; ++v)
        if (alive[v]) {
            long long w = W[v] - lambda;
            if (w >= 0) { add_arc(s, V[v], w); sumPos += w; ++edge_count; }
            else        { add_arc(V[v], t, -w); ++edge_count; }
        }

    // 前驱约束
    for (int i = 0, v = 0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v) {
            if (!alive[v] || V[v] == lemon::INVALID) continue;
            int pre1 = safe_id(i - 1, j);
            int pre2 = safe_id(i, j - 1);
            if (pre1 >= 0 && alive[pre1] && V[pre1] != lemon::INVALID)
                add_arc(V[v], V[pre1], INF), ++edge_count;
            if (pre2 >= 0 && alive[pre2] && V[pre2] != lemon::INVALID)
                add_arc(V[v], V[pre2], INF), ++edge_count;
        }

    if (verbose)
        std::cout << "  graph: nodes=" << (alive_cnt + 2)
                  << " arcs~=" << edge_count << std::endl;

    lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<long long>>
        pf(g, cap, s, t);
    pf.runMinCut(); ++flows_counter;

    std::vector<char> inS(n, 0);
    int sz = 0; long long sumVal = 0;
    for (int v = 0; v < n; ++v)
        if (alive[v] && V[v] != lemon::INVALID) {
            bool onS = pf.minCut(V[v]);
            if (onS) { inS[v] = 1; ++sz; sumVal += (W[v] - lambda); }
        }

    if (verbose)
        std::cout << "  cut: |S|=" << sz
                  << " sumVal=" << sumVal
                  << " sumPos(src)=" << sumPos << std::endl;

    inS_out.swap(inS);
    return {sz, sumVal};
}



// —— 二分 λ 找最大密度闭包（保证非空），返回 inS ——
// predicate(mid): 存在非空闭包 S 使 ∑(W - mid) ≥ 0
// 取最大的 mid 使 predicate 为真；若所有 mid 都假，则取 W 最大的单点闭包。
static std::vector<char>
max_density_closure(int m, int h,
                    const std::vector<char>& alive,
                    const std::vector<long long>& W,
                    long long &flows_counter,
                    bool verbose)
{
    const int n = m*h;
    const long long INF = 1e12;

    int alive_cnt = 0; for (int v=0; v<n; ++v) if (alive[v]) ++alive_cnt;
    if (alive_cnt == 0) return std::vector<char>(n, 0);

    long long lo = std::numeric_limits<long long>::max();
    long long hi = std::numeric_limits<long long>::min();
    for (int v=0; v<n; ++v) if (alive[v]) {
        lo = std::min(lo, W[v]);
        hi = std::max(hi, W[v]);
    }
    // 放宽边界，确保覆盖所有跳变
    lo = lo - (2*alive_cnt + 10);
    hi = hi + (2*alive_cnt + 10);

    std::vector<char> bestS; bestS.reserve(n);

    auto predicate = [&](long long mid)->std::pair<bool, std::vector<char>>{
        std::vector<char> inS;
        auto [sz, val] = max_closure_at_lambda(m,h,alive,W,mid,INF,(int&)flows_counter,inS,verbose);
        bool ok = (sz > 0) && (val >= 0);
        return {ok, std::move(inS)};
    };

    // 二分上中位
    while (lo < hi) {
        long long mid = lo + (hi - lo + 1)/2;
        auto [ok, inS] = predicate(mid);
        if (ok) { lo = mid; bestS = std::move(inS); }
        else    { hi = mid - 1; }
    }

    // 若仍为空（极端全负），退化为取权重 W 最大的单点
    if (bestS.empty()) {
        int argmax = -1; long long bestW = std::numeric_limits<long long>::min();
        for (int v=0; v<n; ++v) if (alive[v]) {
            if (W[v] > bestW) { bestW = W[v]; argmax = v; }
        }
        std::vector<char> S(n, 0);
        if (argmax >= 0) S[argmax] = 1;
        if (verbose) std::cout << "[max_density] fallback single v=" << argmax << " W=" << bestW << std::endl;
        return S;
    }

    // 保证非空闭包
    if (verbose) {
        int sz = 0; for (char c: bestS) if (c) ++sz;
        std::cout << "[max_density] chosen |S|=" << sz << " at lambda=" << lo << std::endl;
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
        vector<char> S = max_density_closure(m,h,A,W,flows,true);


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
        // ===== 递归 Sidney 分解：最大密度闭包放前，再递归补图 =====
    // std::vector<char> alive
    bool verbose = cfg_.verbose;

    // ---- build(...) 递归（稳健版） ----
    std::function<void(const std::vector<char>&)> build =
    [&](const std::vector<char>& A){
        int rest = 0;
for (char c : A) if (c) ++rest;
if (rest <= 0) return;


        if (rest == 1) { // 叶子：唯一一个 alive
            for (int v = 0; v < n; ++v) if (A[v]) { order.push_back(v); break; }
            return;
        }

        if (verbose) {
            std::cout << "[build] rest=" << rest << std::endl;
        }

        // 取最大密度闭包（非空；内部保证非空，即便全负）
        std::vector<char> S = max_density_closure(m,h,A,W,flows,verbose);

        // 统计 S 大小
        int szS = 0; for (char c : S) if (c) ++szS;
        if (szS <= 0) {
            // 极端兜底：若 S 非法为空（理论不应发生），线性取一个点
            if (verbose) std::cout << "[build] S empty fallback\n";
            for (int v=0; v<n; ++v) if (A[v]) { order.push_back(v); return; }
        }

        // 拆分 S 与补
        std::vector<char> A_in(n,0), A_out(n,0);
        for (int v=0; v<n; ++v) if (A[v]) {
            if (S[v]) A_in[v]=1; else A_out[v]=1;
        }

        // 先递归 S（放前），再递归补（放后）
        build(A_in);
        build(A_out);
        if (rest <= 1) {
    for (int v = 0; v < n; ++v)
        if (A[v]) { order.push_back(v); break; }
    return;
}

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
