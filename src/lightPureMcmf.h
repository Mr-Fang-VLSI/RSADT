#pragma once
#include <vector>
#include <cstdint>

struct LightResult {
    int m = 0, h = 0;
    long long dV = 1;
    std::vector<std::vector<long long>> y_order; // ranks 1..m*h
    std::vector<long long> x_of_col;             // all zeros in single-column test
    long long mcmf_cost = 0;
    long long actual_hpwl = 0;
    bool oc_ok = false;
    bool unique_ok = false;
};

class lightPureMcmf {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
    };

    explicit lightPureMcmf(const Config& cfg) : cfg_(cfg) {}
    LightResult solve(int m, int h);

private:
    Config cfg_;

    static long long hpwl_sum(const std::vector<std::vector<long long>>& y,
                              const std::vector<long long>& x,
                              long long dV);
    static bool check_OC(const std::vector<std::vector<long long>>& y);
};
