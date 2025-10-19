#pragma once
#include <vector>
#include <cstdint>

struct KcParamResult {
    int m = 0, h = 0;
    long long dV = 1;

    std::vector<std::vector<int>> y_order; // rank 1..m*h
    std::vector<long long> x_of_col;       // 单列=全0

    long long cost = 0;         // dV * sum w*rank
    long long actual_hpwl = 0;  // 实测 HPWL（四邻边）
    bool oc_ok = false;
    bool unique_ok = false;

    // 统计
    int cuts_total = 0;         // min-cut 次数
};

class KClosureParamOCPlacer {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
        // λ 的整数尺度：所有容量（权与 λ）都用 long long 整数
        long long LAM_SCALE = 1;  // dV 已经是整数时，可取 1
        // 二分最大迭代（权范围很小，20 足够）
        int max_bisect_iter = 24;
    };

    explicit KClosureParamOCPlacer(const Config& cfg): cfg_(cfg) {}
    KcParamResult solve(int m, int h);

private:
    Config cfg_;

    static inline int id(int i, int j, int h) { return i*h + j; }
    static inline int base_weight(int i, int j, int m, int h) {
        int w = 0;
        if (i == 0)     w -= 1;
        if (i == m-1)   w += 1;
        if (j == 0)     w -= 1;
        if (j == h-1)   w += 1;
        return w; // ∈{-2,-1,0,1,2}
    }

    // 一次闭包：给定 λ 与 keep（强制包含），返回源侧 inS
    int max_closure_with_keep(int m, int h,
                              const std::vector<long long>& W, // scaled by dV
                              long long lam,
                              const std::vector<char>& keep,
                              std::vector<char>& inS,
                              int& cuts_counter,
                              bool verbose);

    static long long hpwl_edges_sum(const std::vector<std::vector<int>>& y,
                                    long long dV);
    static bool check_global_OC(const std::vector<std::vector<int>>& y);
};
