#pragma once
#include <vector>
#include <cstdint>
#include <utility>

struct PCResult {
    int m = 0, h = 0;
    long long dV = 1;

    std::vector<std::vector<int>> y_order; // rank 1..m*h
    std::vector<long long> x_of_col;       // 单列=全0

    long long cost = 0;         // 目标值（= dV * sum w*rank）
    long long actual_hpwl = 0;  // 实测 HPWL（四邻边）
    bool oc_ok = false;
    bool unique_ok = false;

    // 统计
    int flows_run = 0;          // 进行的 min-cut 次数
};

class ParamClosureOCPlacer {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
        long long maxT_steps = -1; // 预留：单列 max-T（默认禁用）
    };

    explicit ParamClosureOCPlacer(const Config& cfg): cfg_(cfg) {}
    PCResult solve(int m, int h);

private:
    Config cfg_;

    // 权重：边界 +1/-1，内部 0
    static inline int base_weight(int i, int j, int m, int h) {
        int w = 0;
        if (i == 0)     w -= 1;   // bottom
        if (i == m-1)   w += 1;   // top
        if (j == 0)     w -= 1;   // left
        if (j == h-1)   w += 1;   // right
        return w; // ∈{-2,-1,0,1,2}
    }

    // 真值核验器
    static long long hpwl_edges_sum(const std::vector<std::vector<int>>& y,
                                    const std::vector<long long>& x,
                                    long long dV);
    static bool check_global_OC(const std::vector<std::vector<int>>& y);
};
