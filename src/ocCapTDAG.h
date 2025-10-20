#pragma once
#include <vector>
#include <string>
#include <cstdint>

struct CapTDAGResult {
    int m=0, h=0, n=0, T=0;
    long long dV=1;
    std::vector<std::vector<int>> y_order;
    long long hpwl=0;
    int max_delta=0;
    bool oc_ok=false;
};

class ocCapTDAG {
public:
    struct Config {
        long long dV = 1;         // edge weight scale
        int T = -1;               // max delta, <=0 => no limit (=n)
        int vlevel = 2;           // 1..3
        bool progress = true;     // level progress into log
        std::string logfile = "run.log";
        bool log_append = false;
        int max_labels_per_key = 0; // 0 => no truncation (exact)
    } cfg_;

    explicit ocCapTDAG(const Config& c): cfg_(c) {}

    CapTDAGResult solve(int m, int h);

    // helpers
    static long long hpwl_neighbors(const std::vector<std::vector<int>>& y, long long dV);
    static bool check_OC(const std::vector<std::vector<int>>& y);
    static int  max_delta_adj(const std::vector<std::vector<int>>& y);
    static long long weight_ij(int i,int j,int m,int h);

private:
    void compute_windows_spanT(int m,int h,int T,int n,
                               std::vector<std::vector<int>>& LB,
                               std::vector<std::vector<int>>& UB) const;
};
