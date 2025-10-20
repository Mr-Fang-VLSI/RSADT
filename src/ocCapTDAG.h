#pragma once
#include <vector>
#include <string>
#include <cstdint>

struct CapTDAGResult {
    int m=0, h=0, n=0, T=0;
    long long dV=1;
    std::vector<std::vector<int>> y_order;   // 1..n
    long long hpwl=0;
    bool oc_ok=true;
    int max_delta=0;
};

class ocCapTDAG {
public:
    struct Config {
        long long dV = 1;
        int T = 0;                  // 0 or <=0 => no limit (T=n)
        int vlevel = 2;             // 0~3
        bool progress = true;
        std::string logfile = "run.log";
        bool log_append = false;
        int max_labels_per_key = 64; // 初始 K；层内会自适应增大（≤4096）
    };
    explicit ocCapTDAG(const Config& c): cfg_(c) {}

    CapTDAGResult solve(int m, int h);

private:
    // base-(h+1) encoding
    static inline uint64_t enc_inc(uint64_t key, int idx, const std::vector<uint64_t>& powB) {
        return key + powB[idx];
    }
    static inline int digit_at(uint64_t key, int idx, const std::vector<uint64_t>& powB, uint64_t B){
        return (int)((key / powB[idx]) % B);
    }

    static long long weight_ij(int i,int j,int m,int h);
    static long long hpwl_neighbors(const std::vector<std::vector<int>>& y, long long dV);
    static bool check_OC(const std::vector<std::vector<int>>& y);
    static int  max_delta_adj(const std::vector<std::vector<int>>& y);

    void compute_windows_spanT(int m,int h,int T,int n,
                               std::vector<std::vector<int>>& LB,
                               std::vector<std::vector<int>>& UB) const;

    Config cfg_;
};
