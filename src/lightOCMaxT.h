#pragma once
#include <vector>
#include <string>
#include <cstdint>

struct OCMaxTResult {
    int m=0, h=0, n_full=0, T=0;
    long long dV=1;
    std::vector<std::vector<int>> y_order;
    long long hpwl_full=0;
    long long hpwl_prefix=0;
    long long total_cost=0;
    bool oc_ok=true;
};

class lightOCMaxT {
public:
    struct Config {
        long long dV = 1;
        bool progress = false;        // 是否记录分层进度
        int prefix_strategy = 2;      // 当前使用 2 = Δ≤T 窗口 + Hirschberg
        int omp_threads = 0;          // 0=auto
        int verbose_level = 1;        // 0~3，3 最详细
        std::string logfile = "run.log";
        bool log_append = false;
    };
    explicit lightOCMaxT(const Config& c): cfg_(c) {}

    OCMaxTResult solve(int m, int h, int T_in);

    // Base-(h+1) 编码/解码（m≤32, h≤32）
    static inline uint64_t encode_digit_inc(uint64_t key, int idx, const std::vector<uint64_t>& powB) {
        return key + powB[idx];
    }
    static inline uint64_t encode_digit_dec(uint64_t key, int idx, const std::vector<uint64_t>& powB) {
        return key - powB[idx];
    }
    static inline int digit_at(uint64_t key, int idx, const std::vector<uint64_t>& powB, uint64_t B){
        return (int)((key / powB[idx]) % B);
    }

    // 评估/校验
    static long long hpwl_full_neighbors(const std::vector<std::vector<int>>& y, long long dV);
    static long long hpwl_prefix_neighbors(const std::vector<std::vector<int>>& y, long long dV);
    static bool check_OC_prefix(const std::vector<std::vector<int>>& y);

    Config cfg_;
};
