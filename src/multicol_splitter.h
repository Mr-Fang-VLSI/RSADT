#pragma once
#include <string>
#include <vector>

struct MultiColOptions {
    int m = 0;                  // rows
    int h = 0;                  // cols (logic)
    int K = 2;                  // number of physical DSP columns (典型为 2/4/5)
    int T = 0;                  // target T
    int rounds = 100;           // outer iterations
    std::string policy = "soft"; // "soft" | "lagrange" | "irl1"
    std::string single_bin = "build/oc_policy_singlecol"; // underlying solver binary

    // soft policy
    double alpha = 0.9, eta = 1.0, gamma = 0.3, beta = 3.0;
    int homotopy = 1;           // T' tightening on/off
    double sigma = 0.8, shrink = 0.9;

    // lagrange
    double rho = 0.3;

    // irl1
    double lam = 2.0, delta = 2.0;

    // shared backend options
    int freeze = 1;
    int dV = 1;
    int pow2 = 1;
    int progress = 1;
    long long cap = 64;
    long long scale = 1000;

    // NEW: snake (logical left-right flip) on odd strips.
    // Implementation: for odd strip s, index Y as y_base[i][w-1-j_local] when emitting output.
    // This is a pure "parser-time" trick; we do not change DP or weights.
    int snake = 0;              // 0=off (default), 1=flip logical j on odd strips

    // working dir
    std::string work_dir = ".";
};

class MultiColSplitter {
public:
    // Run single-column solver on one strip (m x (h/K)) and replicate to K columns.
    // Returns path to final placement file (placement{m}_Col_{K}_T_{T}.txt) on success, empty string otherwise.
    std::string run(const MultiColOptions& opt);

private:
    static bool ensure_dir(const std::string& dir);
    static int  run_cmd(const std::string& cmd);

    // Parse "placement{m}_Col_1_T_{T}.txt" into y_base[i][j_local] (only Y used; X is 0 in single-col).
    static bool parse_singlecol_placement(const std::string& path,
                                          int m, int w,
                                          std::vector<std::vector<long long>>& y_base);

    // Write final multi-col placement: replicate base Y to each strip, X = strip id
    // If snake=1, for odd strips s, use j_idx = w-1-j_local when indexing y_base (logical left-right flip).
    static bool write_multicol_placement(const std::string& out_path,
                                         int m, int h, int K, int w,
                                         const std::vector<std::vector<long long>>& y_base,
                                         int snake);
};
