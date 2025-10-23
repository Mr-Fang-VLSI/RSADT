#pragma once
#include <vector>
#include <cstdint>
#include <string>

struct DpStats {
    double time_sec = 0.0;
    double peak_rss_mb = 0.0;
    size_t final_states = 0;
};

struct DpResult {
    long long best_hpwl = 0;
    std::vector<std::vector<int>> y;  // ?? reconstruct=true ???
    DpStats stats;
};

class DpOcSolver {
public:
    struct Config {
        bool reconstruct = false;   // ???? y
        bool verbose = false;       // ????
    };
    explicit DpOcSolver(const Config& cfg): cfg_(cfg) {}

    // ???:?? Max-T?????? OC + AllDifferent ?? HPWL
    DpResult solve(int m, int h);

    // ??/??
    static long long hpwl_neighbors(const std::vector<std::vector<int>>& y);
    static long long hpwl_boundary(const std::vector<std::vector<int>>& y);

private:
    Config cfg_;

    // ?????(wV=wH=1)
    static void build_C_unit(int m,int h, std::vector<std::vector<long long>>& C);

    // Base-(h+1) ????
    static std::vector<uint64_t> powB_list(int m, uint64_t B);
    static inline uint64_t encode_inc(uint64_t key, int idx, const std::vector<uint64_t>& pw){
        return key + pw[idx];
    }
    static inline int digit_at(uint64_t key, int i, const std::vector<uint64_t>& pw, uint64_t B){
        return int((key / pw[i]) % B);
    }

    // ????(ru_maxrss)
    static double get_peak_rss_mb();
};
