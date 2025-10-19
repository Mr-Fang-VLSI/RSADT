#pragma once
#include <vector>
#include <cstdint>
#include <utility>
#include <string>

/**
 * IDA*-based optimal single-column placer under OC (order consistency).
 * - Grid: m rows (0..m-1, bottom to top), h columns (0..h-1, left to right).
 * - We assign each cell (i,j) a distinct rank 1..m*h such that if i1<=i2 and j1<=j2 then
 *   rank(i1,j1) <= rank(i2,j2).  (Global OC.)
 * - Objective: minimize the HPWL sum over all horizontal and vertical neighbor nets
 *   when all columns are mapped to a single physical site column (x=0).
 *   With OC, the sum telescopes and equals sum_{i,j} w(i,j) * rank(i,j) * dV
 *   where w(i,j) = (top?+1 : bottom?-1 : 0) + (right?+1 : left?-1 : 0).
 *   Corners become: TR:+2, BL:-2, TL:0, BR:0; internal cells:0; edges: +/-1.
 *
 * This solver guarantees a global optimum (no pruning heuristics that risk suboptimality).
 * It uses IDA* with an admissible lower bound based on the multiset of remaining weights.
 */
struct IdaResult {
    int m = 0, h = 0;
    long long dV = 1;
    // y_order[i][j] in {1..m*h}; physical y = y_order * dV
    std::vector<std::vector<long long>> y_order;
    // all columns mapped to a single physical site column ⇒ x_of_col = 0
    std::vector<long long> x_of_col;

    long long objective_cost = 0;   // sum w*rank*dV
    long long actual_hpwl = 0;      // recomputed HPWL from y_order (sanity check)
    bool oc_ok = false;             // global OC check
    bool unique_ok = false;         // all ranks unique
    // search stats
    uint64_t nodes_expanded = 0;
    uint64_t dfs_calls = 0;
    uint64_t ida_iterations = 0;
};

class IdaStarColumnPlacer {
public:
    struct Config {
        long long dV = 1;      // vertical unit distance
        bool verbose = true;   // print progress
        // Optional hard limits; 0 means unlimited within process constraints.
        uint64_t max_nodes = 0;
        uint64_t max_dfs_calls = 0;
    };
    explicit IdaStarColumnPlacer(const Config& cfg) : cfg_(cfg) {}

    // Solve single-column placement (all logical columns map to one site column).
    IdaResult solve(int m, int h);

private:
    Config cfg_;

    // === utilities ===
    static inline int weight_of(int i, int j, int m, int h) {
        int w = 0;
        if (i == m-1) w += 1; else if (i == 0) w -= 1;
        if (j == h-1) w += 1; else if (j == 0) w -= 1;
        return w; // in {-2,-1,0,1,2}
    }
    static long long hpwl_edges_sum(const std::vector<std::vector<long long>>& y,
                                    const std::vector<long long>& x,
                                    long long dV);
    static bool check_global_OC(const std::vector<std::vector<long long>>& y);
    static bool check_unique_1_to_N(const std::vector<std::vector<long long>>& y);

    // === IDA* core ===
    struct State {
        std::vector<int> a;     // frontier counts per row (0..h), non-increasing: a[m-1] <= ... <= a[0]
        int placed = 0;         // how many cells placed (also next rank-1)
        long long g = 0;        // accumulated cost
        // counts of remaining weights w in order [-2,-1,0,1,2] -> idx 0..4
        int rem_cnt[5] = {0,0,0,0,0};
    };

    long long initial_bound(const State& s, int m, int h) const;
    long long heuristic_lb(const State& s, int m, int h) const;
    bool ida_search(long long& bound, State& s,
                    const int m, const int h,
                    std::vector<std::pair<int,int>>& picks,
                    IdaResult& result);

    bool dfs(long long bound, long long& next_bound, State& s,
             const int m, const int h,
             std::vector<std::pair<int,int>>& picks,
             IdaResult& result);
};
