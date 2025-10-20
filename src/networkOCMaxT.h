#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include "lightOCMaxT.h"
// 与 lightOCMaxT.h 中的结构保持一致（用于统一输出结果）
// struct OCMaxTResult {
//     int m=0, h=0, n_full=0, T=0;
//     long long dV=1;
//     std::vector<std::vector<int>> y_order;
//     long long hpwl_full=0;
//     long long hpwl_prefix=0;
//     long long total_cost=0;
//     bool oc_ok=true;
// };

// ============ RCDC 完整 DAG 版本 ============
// 说明：该类实现“网络化”的完整层序最短路。
// 状态包含：天际线 a、行时钟 R[i]（该行最近一次放置层号）
// 和列时钟 C[j]（该列竖直父的实际层号），从而在每次转移时
// 精确检查 Δ≤T（避免 post-check 失败）。
class networkOCMaxT {
public:
    struct Config {
        long long dV = 1;           // 邻接边权的缩放（等价 HPWL 的系数）
        bool progress = true;       // 是否打印分层统计到日志
        int  verbose_level = 2;     // 0~3
        std::string logfile = "run.log";
        bool log_append = false;    // true: 追加; false: 覆盖
    };
    explicit networkOCMaxT(const Config& c): cfg_(c) {}

    // 求解入口：m 行、h 列、Δ≤T
    OCMaxTResult solve(int m, int h, int T_in);

private:
    Config cfg_;
};
