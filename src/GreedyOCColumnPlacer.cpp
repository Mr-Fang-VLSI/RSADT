#include "GreedyOCColumnPlacer.h"
#include <queue>
#include <unordered_set>
#include <iostream>
#include <algorithm>

using std::vector;
using std::cout;
using std::endl;

long long GreedyOCColumnPlacer::hpwl_edges_sum(const vector<vector<int>>& y,
                                               const vector<long long>& x,
                                               long long dV) {
    (void)x; // 单列 x 恒为 0
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long S = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j+1 < h; ++j)
            S += (long long)std::llabs(y[i][j+1]-y[i][j]) * dV; // 水平邻边（dx=0）
    for (int i = 0; i+1 < m; ++i)
        for (int j = 0; j < h; ++j)
            S += (long long)std::llabs(y[i+1][j]-y[i][j]) * dV; // 竖向邻边
    return S;
}

bool GreedyOCColumnPlacer::check_global_OC(const vector<vector<int>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
    for (int j1 = 0; j1 < h; ++j1)
    for (int i2 = i1; i2 < m; ++i2)
    for (int j2 = j1; j2 < h; ++j2)
        if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

GreedyResult GreedyOCColumnPlacer::solve(int m, int h) {
    const int N = m * h;
    GreedyResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order.assign(m, vector<int>(h, 0));
    R.x_of_col.assign(h, 0);

    // 前沿: a[i] 表示第 i 行已经放置的列数（0..h），保持非增
    vector<int> a(m, 0);

    // 候选最大堆（按 w 降序；若 w 相等，次序用 (i,j) 的一个稳定规则，尽量兼容你之前对称偏好）
    struct Cand {
        int w, i, j; // j = a[i]
        int t;       // 选择该候选时的 rank（只用于打印/调试；排序以 w 为主）
    };
    struct Cmp {
        bool operator()(const Cand& A, const Cand& B) const {
            if (A.w != B.w) return A.w < B.w;  // 最大堆：w 大者先
            // 次关键字（稳定一些）：优先低 j（左），再低 i（下）
            if (A.j != B.j) return A.j > B.j;
            return A.i > B.i;
        }
    };
    std::priority_queue<Cand, std::vector<Cand>, Cmp> pq;

    auto try_push = [&](int i, int t_for_dbg){
        if (i < 0 || i >= m) return;
        if (!eligible(i, a, h)) return;
        int j = a[i];
        int w = weight(i, j, m, h);
        pq.push(Cand{w, i, j, t_for_dbg});
    };

    // 初始化：所有可扩展行各推下一个格
    for (int i = 0; i < m; ++i) try_push(i, 0);

    long long Z = 0;
    int placed = 0;

    while (placed < N) {
        if (pq.empty()) {
            throw std::runtime_error("GreedyOCColumnPlacer: empty frontier (should not happen).");
        }
        Cand c = pq.top(); pq.pop();

        // 懒惰校验：候选必须与当前前沿一致，否则丢弃
        if (!(c.j == a[c.i] && eligible(c.i, a, h))) continue;

        int i = c.i, j = c.j;
        int t = placed + 1;

        // 赋 rank
        R.y_order[i][j] = t;
        // 累计目标值（等价 HPWL 的线性形式）
        Z += (long long)c.w * (long long)t * cfg_.dV;

        // 更新前沿
        a[i] += 1;
        placed += 1;

        // 该行下一候选
        try_push(i, t);
        // 可能解锁下一行
        try_push(i+1, t);

        // TODO: 将来支持 maxT: 需要对 (i,j-1)->(i,j) 的 “间隔计数” 做上界检查
        // 若某行 (i) 的相邻列差可能超出 T，在入堆时把“必须选”的候选放到专门的集合（或抬高 w 到 +∞）
    }

    R.cost = Z;
    R.actual_hpwl = hpwl_edges_sum(R.y_order, R.x_of_col, cfg_.dV);
    R.oc_ok = check_global_OC(R.y_order);

    // 唯一占用 (x,y)
    {
        std::unordered_set<uint64_t> occ; occ.reserve((size_t)N*2);
        bool ok = true;
        for (int ii = 0; ii < m; ++ii)
            for (int jj = 0; jj < h; ++jj) {
                uint64_t x = 0ull;
                uint64_t yv = (uint64_t)R.y_order[ii][jj];
                uint64_t key64 = (x << 32) ^ yv;
                if (!occ.insert(key64).second) ok = false;
            }
        R.unique_ok = ok;
    }
    R.steps = placed;

    if (cfg_.verbose) {
        cout << "[Greedy-OC] placed=" << R.steps
             << ", cost=" << R.cost
             << ", HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok ? "OK":"FAIL")
             << ", Unique=" << (R.unique_ok ? "OK":"FAIL") << endl;
    }
    return R;
}
