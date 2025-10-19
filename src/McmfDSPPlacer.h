#pragma once
#include <vector>
#include <string>
#include <utility>
#include <limits>
#include <cstdint>

struct Segment {
    int start_col = 0;
    int width = 0;
    int g = 1;
    bool mirrored = false; // not used when enforce_global_oc=true
};

struct PlacementResult {
    int m = 0, h = 0;
    int used_site_cols = 0;
    long long dH = 10, dV = 1;

    // Physical y used for HPWL (after minimal baseline shifts to satisfy OC)
    std::vector<std::vector<long long>> y_order; // size m x h
    // Global display labels 1..(m*h) (for debug only)
    std::vector<std::vector<long long>> y_label; // size m x h

    std::vector<long long> x_of_col;

    std::vector<Segment> segments;

    long long mcmf_cost = 0;
    long long actual_hpwl = 0;

    bool oc_ok = false;
    bool unique_ok = false;

    // Max-T report
    long long max_step_observed = 0; // in steps (before multiplying dV)
    long long max_T_limit = -1;      // in steps; <0 means disabled
    bool maxT_ok = true;
};

class McmfDSPPlacer {
public:
    struct Config {
        long long dH = 10;
        long long dV = 1;
        int max_site_cols = 100;
        bool restrict_equal_width = true;
        bool verbose = true;
        bool enforce_global_oc = true;   // global nondecreasing y across columns
        long long max_T_steps = -1;      // <=0: disabled; else: max per-horizontal-edge vertical steps
    };

    explicit McmfDSPPlacer(const Config& cfg);
    PlacementResult solve(int m, int h);

private:
    Config cfg_;

    // ---- RSAD cost & mapping ----
    std::pair<long long,int> LI_closed_form_int(int m, int w) const; // min over g
    long long LI_by_g_int(int m, int w, int g) const;                // evaluate at g
    std::vector<std::vector<long long>> rsad_mapping(int m, int w, int g) const;

    // Vertical penalty between two patterns (equal width) used in cost
    long long vpen_between_patterns(int m, int w, int g_prev, int g_next) const;

    // ---- max-T helpers ----
    long long max_intra_step(int m, int w, int g) const;                          // max |y(i,j+1)-y(i,j)|
    long long max_inter_step_after_shift(int m, int w, int g_prev, int g_next) const; // = range(R-L)

    // ---- MCMF (OC-aware) ----
    // Choose equal width w and a g for each segment; include Vpen on transitions and filter by max-T.
    std::tuple<int,int,long long,std::vector<int>> solve_mcmf_equal_width_oc(int m, int h) const;

    // ---- Utilities ----
    static long long hpwl_edges_sum(const std::vector<std::vector<long long>>& y,
                                    const std::vector<long long>& x,
                                    long long dV);
    static bool check_global_OC(const std::vector<std::vector<long long>>& y);
    static long long max_horizontal_step_overall(const std::vector<std::vector<long long>>& y);

    // RSAD region helpers
    static long long OLL(int i, int j);
    static long long OLR(int i, int j, int g);
    static long long OUL(int i, int j, int g);
    static long long OUR(int i, int j, int g);
};
