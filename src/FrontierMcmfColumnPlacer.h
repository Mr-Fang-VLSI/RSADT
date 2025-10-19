#pragma once
#include <vector>
#include <cstdint>
#include <string>


struct FrontierResult {
    int m = 0, h = 0;
    long long dV = 1;

    // 全局排名（1..m*h），满足 OC；物理 y = rank * dV；物理 x = 0（单列）
    std::vector<std::vector<int>> y_order;
    std::vector<long long> x_of_col; // 全 0（单列）

    long long mcmf_cost = 0;    // 与 HPWL 精确相等（单位：物理长度）
    long long actual_hpwl = 0;

    bool oc_ok = false;
    bool unique_ok = false;

    // 调试统计
    size_t states_expanded = 0;
    size_t states_cached = 0;
};

class FrontierMcmfColumnPlacer {
public:
    struct Config {
        long long dV = 1;   // 竖向单位长度
        bool verbose = true;
    };

    explicit FrontierMcmfColumnPlacer(const Config& cfg): cfg_(cfg) {}

    // 求解（单列放置：dH→∞ 场景）
    FrontierResult solve(int m, int h);

private:
    Config cfg_;

    // === 代价权重：边界望远镜分解（与 HPWL 等价） ===
    static inline int weight(int i, int j, int m, int h) {
        int w = 0;
        if (i == 0)     w -= 1;   // bottom
        if (i == m-1)   w += 1;   // top
        if (j == 0)     w -= 1;   // left
        if (j == h-1)   w += 1;   // right
        return w; // ∈{-2,-1,0,1,2}
    }

    // === 真值核验器 ===
    static long long hpwl_edges_sum(const std::vector<std::vector<int>>& y,
                                    const std::vector<long long>& x,
                                    long long dV);
    static bool check_global_OC(const std::vector<std::vector<int>>& y);
};
