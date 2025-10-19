#pragma once
#include <vector>
#include <cstdint>

struct OCShortestResult {
    int m = 0, h = 0, n = 0;
    long long dV = 1;
    // y_order[i][j] = rank 1..n
    std::vector<std::vector<int>> y_order;
    long long total_cost = 0;
    long long hpwl = 0;
    bool oc_ok = false;
};

class lightOCShortest {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
    };
    explicit lightOCShortest(const Config& cfg) : cfg_(cfg) {}
    OCShortestResult solve(int m, int h);

private:
    Config cfg_;

    static long long hpwl_sum(const std::vector<std::vector<int>>& y, long long dV);
    static bool check_OC(const std::vector<std::vector<int>>& y);
};
