#pragma once
#include <vector>
#include <cstdint>

struct OCMaxTResult {
    int m = 0, h = 0, n_full = 0;   // n_full = m*h
    int T = 0;                      // 实际使用的层数（钳制后的）
    long long dV = 1;

    // y_order[i][j]：若 1..T 内被赋值则为其秩，否则为 0（当 T<n_full）
    std::vector<std::vector<int>> y_order;

    // 目标值（分层 DAG 0..T 的最短路成本）
    long long total_cost = 0;

    // 仅当 T==n_full 时：完整 HPWL（与 total_cost 一致）
    long long hpwl_full = 0;

    // 当 T<n_full 时：只统计两端都已赋秩的相邻边
    long long hpwl_prefix = 0;

    bool oc_ok = false;
};

class lightOCMaxT {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
        bool progress = false;
        int  omp_threads = 0;   // 0 -> all cores
        bool low_mem = true;    // 使用 Hirschberg 低内存回溯
    };
    explicit lightOCMaxT(const Config& cfg) : cfg_(cfg) {}

    // 主入口：m,h, 以及 Max-T（可>min(m,h)，最终 T 会被钳制到 [1, m*h]）
    OCMaxTResult solve(int m, int h, int T);

    // === 工具函数（类外也要用，公开） ===
    static inline uint64_t encode_digit_inc(uint64_t key, int idx,
                                            const std::vector<uint64_t>& powB) {
        return key + powB[idx];
    }
    static inline uint64_t encode_digit_dec(uint64_t key, int idx,
                                            const std::vector<uint64_t>& powB) {
        return key - powB[idx];
    }
    static inline int digit_at(uint64_t key, int idx,
                               const std::vector<uint64_t>& powB, uint64_t B){
        return int((key / powB[idx]) % B);
    }

private:
    Config cfg_;

    // 邻接 HPWL（完整网格）
    static long long hpwl_full_neighbors(const std::vector<std::vector<int>>& y, long long dV);
    // 前缀 HPWL：仅统计两端都>0 的相邻边
    static long long hpwl_prefix_neighbors(const std::vector<std::vector<int>>& y, long long dV);
    static bool check_OC(const std::vector<std::vector<int>>& y);
};
