#pragma once
#include "lightOCMaxT.h"   // 复用 OCMaxTResult 定义
#include <string>

class networkOCMaxT {
public:
    struct Config {
        long long dV = 1;            // HPWL 等价的系数（邻接边权缩放）
        bool progress = true;         // 输出层级统计
        int  verbose_level = 2;       // 0~3
        std::string logfile = "run.log";
        bool log_append = false;

        // 性能/稳健选项
        int  log_every_levels = 1;    // v>=3 时每多少层打印一次（1=每层）
        int  cap_per_bucket   = 0;    // 每个 akey 的非支配前沿上限（0=不开）
        bool enable_dominance = true; // 打开/关闭主导裁剪
    };
    explicit networkOCMaxT(const Config& c): cfg_(c) {}

    OCMaxTResult solve(int m, int h, int T_in);

private:
    Config cfg_;
};
