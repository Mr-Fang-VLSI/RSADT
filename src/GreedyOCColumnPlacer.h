#pragma once
#include <vector>
#include <cstdint>

struct GreedyResult {
    int m = 0, h = 0;
    long long dV = 1;

    // 全局排名: 1..(m*h)，满足 OC；物理 y = rank * dV；物理 x 全 0（单列）
    std::vector<std::vector<int>> y_order;
    std::vector<long long> x_of_col;

    long long cost = 0;            // 贪心累计的目标值（单位：物理长度）
    long long actual_hpwl = 0;     // 实算 HPWL（四邻边）
    bool oc_ok = false;
    bool unique_ok = false;

    // 统计
    int steps = 0;
};

class GreedyOCColumnPlacer {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
        // 预留：将来可加入 maxT_steps（单列时：限制 y(i,j+1)-y(i,j) ≤ T）
        long long maxT_steps = -1; // <0 表示不启用
    };

    explicit GreedyOCColumnPlacer(const Config& cfg): cfg_(cfg) {}

    GreedyResult solve(int m, int h);

private:
    Config cfg_;

    static inline int weight(int i, int j, int m, int h) {
        int w = 0;
        if (i == 0)     w -= 1;   // bottom
        if (i == m-1)   w += 1;   // top
        if (j == 0)     w -= 1;   // left
        if (j == h-1)   w += 1;   // right
        return w; // ∈ {-2,-1,0,1,2}
    }

    static bool eligible(int i, const std::vector<int>& a, int h) {
        if (a[i] >= h) return false;
        if (i == 0)    return true;
        return a[i] + 1 <= a[i-1];
    }

    static long long hpwl_edges_sum(const std::vector<std::vector<int>>& y,
                                    const std::vector<long long>& x,
                                    long long dV);
    static bool check_global_OC(const std::vector<std::vector<int>>& y);
};
