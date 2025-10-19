#pragma once
#include <vector>
#include <cstdint>

struct OCShortestResult {
    int m = 0, h = 0, n = 0;
    long long dV = 1;
    std::vector<std::vector<int>> y_order; // rank in [1..n]
    long long total_cost = 0;
    long long hpwl = 0;
    bool oc_ok = false;
};

class lightOCShortest {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
        int  max_steps = -1;   // -1 -> run full n = m*h
        bool progress  = false;
        int  omp_threads = 0;  // 0 -> all cores
        bool low_mem = true;   // ✅ 启用 Hirschberg 低内存回溯
    };
    explicit lightOCShortest(const Config& cfg) : cfg_(cfg) {}
    OCShortestResult solve(int m, int h);
    static inline uint64_t encode_digit_inc(uint64_t key, int idx, const std::vector<uint64_t>& powB) {
        return key + powB[idx];
    }
    static inline uint64_t encode_digit_dec(uint64_t key, int idx, const std::vector<uint64_t>& powB) {
        return key - powB[idx];
    }
    static inline int digit_at(uint64_t key, int idx, const std::vector<uint64_t>& powB, uint64_t B){
        return int((key / powB[idx]) % B);
    }
private:
    Config cfg_;
    static inline long long llabsll(long long x){ return x>=0?x:-x; }
    static long long hpwl_sum(const std::vector<std::vector<int>>& y, long long dV);
    static bool check_OC(const std::vector<std::vector<int>>& y);

    // base-(h+1) 编码的增量与取位
    
};
