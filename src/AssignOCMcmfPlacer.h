#pragma once
#include <vector>
#include <cstdint>

struct AssignOCResult {
    int m = 0, h = 0;
    long long dV = 1;

    std::vector<std::vector<int>> y_order; // rank 1..m*h
    std::vector<long long> x_of_col;       // 单列=全0

    long long cost = 0;         // dV * sum w*rank
    long long actual_hpwl = 0;  // 实测 HPWL（四邻边）
    bool oc_ok = false;
    bool unique_ok = false;

    // 统计
    int nodes = 0;
    long long arcs = 0;
};

class AssignOCMcmfPlacer {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
    };

    explicit AssignOCMcmfPlacer(const Config& cfg): cfg_(cfg) {}
    AssignOCResult solve(int m, int h);

private:
    Config cfg_;

    static inline int vid(int i, int j, int h) { return i*h + j; }
    static inline int base_weight(int i, int j, int m, int h) {
        int w = 0;
        if (i == 0)     w -= 1;
        if (i == m-1)   w += 1;
        if (j == 0)     w -= 1;
        if (j == h-1)   w += 1;
        return w; // ∈{-2,-1,0,1,2}
    }

    static long long hpwl_edges_sum(const std::vector<std::vector<int>>& y,
                                    long long dV);
    static bool check_global_OC(const std::vector<std::vector<int>>& y);
};
