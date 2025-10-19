#pragma once
#include <vector>
#include <cstdint>

struct KClosureResult {
    int m = 0, h = 0;
    long long dV = 1;
    std::vector<std::vector<int>> y_order; // 1..m*h
    long long cost = 0;         // dV * sum w*rank
    long long actual_hpwl = 0;  // 实测 HPWL（四邻边）
    bool oc_ok = false;
    bool unique_ok = false;
    // 统计
    long long mincut_calls = 0;
};

class KClosureExactPlacer {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
    };
    explicit KClosureExactPlacer(const Config& cfg): cfg_(cfg) {}
    KClosureResult solve(int m, int h);

private:
    Config cfg_;
    static inline int vid(int i, int j, int h) { return i*h + j; }
    static inline int wbase(int i, int j, int m, int h) {
        int w = 0;
        if (i==0)   w -= 1;
        if (i==m-1) w += 1;
        if (j==0)   w -= 1;
        if (j==h-1) w += 1;
        return w;
    }
    static long long hpwl_sum(const std::vector<std::vector<int>>& y, long long dV);
    static bool check_OC(const std::vector<std::vector<int>>& y);

    // 一次最大闭包：给定 λ 与 Keep，返回 S（bool 标记）
    std::vector<char> max_closure_lambda(
        int m, int h,
        const std::vector<char>& keep,               // 必须包含的点（上一层 S 的锁定）
        long long lambda,
        long long& cut_pos_sum                       // 返回正边容量和（调试用）
    );

    // 在“0‑边际诱导子图”里，求恰好补 x 个的闭包（保持 keep），精确补齐到 |S|=target
    std::vector<char> closure_fill_zero_margin(
        int m, int h,
        const std::vector<char>& baseS,              // 先验 S（含正边与部分0边）
        const std::vector<char>& zero_mask,          // p(v,λ)==0 的可候选集合
        int target                                    // 期望最终 |S|
    );
};
