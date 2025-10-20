#pragma once
#include "lightOCMaxT.h"   // 直接复用 OCMaxTResult 定义
#include <string>

class networkOCMaxT {
public:
    struct Config {
        long long dV = 1;            // 邻接边权缩放（HPWL 等价因子）
        bool progress = true;         // 是否输出层级统计
        int  verbose_level = 2;       // 0~3；3 打每层
        std::string logfile = "run.log";
        bool log_append = false;

        // 性能/稳健选项
        int  log_every_levels = 1;    // v>=3 时每几层记一条（1=每层）
        int  cap_per_bucket   = 0;    // 每个 akey 的非支配前沿上限（0=不开）
        bool enable_dominance = true; // 主导裁剪开关
    };
    explicit networkOCMaxT(const Config& c): cfg_(c) {}

    OCMaxTResult solve(int m, int h, int T_in);

private:
    Config cfg_;
};
