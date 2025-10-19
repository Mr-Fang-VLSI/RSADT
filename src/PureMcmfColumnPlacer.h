#pragma once
#include <vector>
#include <string>
#include <utility>
#include <cstdint>

struct PureResult {
    int m = 0, h = 0;
    long long dV = 1;

    // y_order[i][j] = 1..(m*h)，满足全局 OC；物理 y = y_order * dV
    std::vector<std::vector<long long>> y_order;

    // 物理 x 坐标（本验证版把所有列打到同一根 site 列 ⇒ 全 0）
    std::vector<long long> x_of_col;

    long long mcmf_cost = 0;     // 来自 LEMON 的 totalCost（单位：物理长度）
    long long actual_hpwl = 0;   // 枚举相邻边实算 HPWL
    bool oc_ok = false;          // 全局 OC 检查
    bool unique_ok = false;      // (x,y) 唯一占用
};

class PureMcmfColumnPlacer {
public:
    struct Config {
        long long dV = 1;      // 竖向单位长度
        bool verbose = true;   // 打印中间信息
    };

    explicit PureMcmfColumnPlacer(const Config& cfg) : cfg_(cfg) {}
    // 仅做“单列”放置（所有逻辑列映射到同一根 site 列），用于 dH→∞ 的验证口径
    PureResult solve(int m, int h);

private:
    Config cfg_;

    // === 工具 ===
    static long long hpwl_edges_sum(const std::vector<std::vector<long long>>& y,
                                    const std::vector<long long>& x,
                                    long long dV);
    static bool check_global_OC(const std::vector<std::vector<long long>>& y);
};
