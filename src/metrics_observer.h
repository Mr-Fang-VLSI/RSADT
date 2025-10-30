#pragma once
#include <vector>
#include <string>
#include <utility>
#include <ostream>

template<typename T> using Mat = std::vector<std::vector<T>>;

struct WeightStats {
    double mean = 0.0, stddev = 0.0, cv = 0.0, median = 0.0, gini = 0.0, top_share = 0.0;
    long long minv = 0, maxv = 0;
    size_t n = 0;
};

struct EdgeLoc {
    bool horiz = true; // true: H(i,j)-(i,j+1); false: V(i,j)-(i+1,j)
    int  i = 0, j = 0;
    double cx = 0.0, cy = 0.0; // center for jump metric
    long long L = 0;
    bool valid = false;
};

struct EdgeStats {
    long long maxL = 0;
    EdgeLoc   maxEdge;
    long long viol_count = 0;
    long long viol_sum = 0; // equals TNS if computed with same T
};

struct JumpStats {
    double manhattan = 0.0;
    double euclid = 0.0;
    bool   orient_flip = false;
};

namespace metrics {

// weights: compute stats over concatenated H and V integer weights
WeightStats compute_weight_stats(const Mat<long long>& WH, const Mat<long long>& WV, double top_ratio = 0.10);

// scan y (1-based) to compute maxL / max-edge location & violations
EdgeStats   compute_edge_stats(const Mat<int>& y, int m, int h, long long T);

// jump distance (curr vs prev)
JumpStats   jump_between(const EdgeLoc& prev, const EdgeLoc& curr, bool has_prev);

// logging helpers
void        write_header(std::ostream& os);
void        write_round(std::ostream& os, int round, long long T,
                        long long hpw_eq, long long hpw_w,
                        const EdgeStats& es, long long WNS, long long TNS, int picked,
                        const WeightStats& ws, const JumpStats& js,
                        double ms_ms, bool feasible, double avg_jump);

} // namespace metrics
