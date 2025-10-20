#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <unordered_map>
#include "lightOCMaxT.h"  // 复用 OCMaxTResult 与配置风格

class networkOCMaxT {
public:
    struct Config {
        long long dV = 1;
        bool progress = true;
        int verbose_level = 2;      // 0~3
        std::string logfile = "run.log";
        bool log_append = true;
        int max_levels_log_step = 16; // 每隔多少层打一次层宽
        int T_shift_guard = 0;      // 预留（本实现不需要 shift）
    };
    explicit networkOCMaxT(const Config& c): cfg_(c) {}

    // 主入口
    OCMaxTResult solve(int m, int h, int T_in);

private:
    // ========== 编码/解码 ==========
    static inline uint64_t encode_inc(uint64_t key, int idx, const std::vector<uint64_t>& powB) {
        return key + powB[idx];
    }
    static inline int digit_at(uint64_t key, int idx, const std::vector<uint64_t>& powB, uint64_t B){
        return (int)((key / powB[idx]) % B);
    }

    // ========== 评估 ==========
    static inline long long weight_ij(int i,int j,int m,int h){
        long long w=0;
        if(i==0)   w -= 1;
        if(i==m-1) w += 1;
        if(j==0)   w -= 1;
        if(j==h-1) w += 1;
        return w;
    }

    // ========== 窗口闭包（与现有实现等价） ==========
    static void compute_windows_spanT(
        int m,int h,int T, int n,
        std::vector<std::vector<int>>& LB,
        std::vector<std::vector<int>>& UB,
        int vlevel,
        const std::string& logpath,
        bool append);

    Config cfg_;
};
