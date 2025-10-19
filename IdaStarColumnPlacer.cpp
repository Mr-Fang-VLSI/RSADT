#include "IdaStarColumnPlacer.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <unordered_set>

using std::vector;
using std::pair;
using std::cout;
using std::endl;

static inline long long sum_arith(long long L, long long R) {
    // sum of integers from L to R inclusive; if L>R returns 0
    if (L > R) return 0;
    return (R - L + 1) * (L + R) / 2;
}

long long IdaStarColumnPlacer::hpwl_edges_sum(const vector<vector<long long>>& y,
                                              const vector<long long>& /*x*/,
                                              long long dV) {
    const int m = (int)y.size();
    const int h = (int)y[0].size();
    long long S = 0;
    // horizontal edges (rows)
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j+1 < h; ++j) {
            long long dy = std::llabs(y[i][j+1] - y[i][j]) * dV;
            S += dy;
        }
    }
    // vertical edges (columns)
    for (int i = 0; i+1 < m; ++i) {
        for (int j = 0; j < h; ++j) {
            long long dy = std::llabs(y[i+1][j] - y[i][j]) * dV;
            S += dy;
        }
    }
    return S;
}

// Global OC: if i1<=i2 and j1<=j2 then y[i1][j1] <= y[i2][j2]
bool IdaStarColumnPlacer::check_global_OC(const vector<vector<long long>>& y) {
    const int m = (int)y.size();
    const int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
        for (int j1 = 0; j1 < h; ++j1)
            for (int i2 = i1; i2 < m; ++i2)
                for (int j2 = j1; j2 < h; ++j2)
                    if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

bool IdaStarColumnPlacer::check_unique_1_to_N(const vector<vector<long long>>& y) {
    const int m = (int)y.size();
    const int h = (int)y[0].size();
    const int N = m*h;
        std::unordered_set<long long> seen;
    seen.reserve((size_t)N*2);
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            long long v = y[i][j];
            if (v < 1 || v > N) return false;
            if (!seen.insert(v).second) return false;
        }
    return true;
}

// Lower bound from multiset of remaining weights ignoring OC.
// We sort weights in descending order and assign to ranks (placed+1..N) ascending.
// This is the rearrangement inequality; gives admissible (optimistic) LB.
long long IdaStarColumnPlacer::heuristic_lb(const State& s, int /*m*/, int /*h*/) const {
        // But we don't need N via that hack. Compute directly:
    int placed = s.placed;
    int total_remaining = s.rem_cnt[0] + s.rem_cnt[1] + s.rem_cnt[2] + s.rem_cnt[3] + s.rem_cnt[4];
    int Ntot = placed + total_remaining;

    // Map weight to count: idx 0..4 -> w = -2,-1,0,1,2
    int c_n2 = s.rem_cnt[0];
    int c_n1 = s.rem_cnt[1];
    int c_0  = s.rem_cnt[2];
    int c_p1 = s.rem_cnt[3];
    int c_p2 = s.rem_cnt[4];

    long long lb = 0;
    long long L = placed + 1;

    // positives at the earliest ranks
    if (c_p2 > 0) {
        long long R = L + c_p2 - 1;
        lb += 2LL * sum_arith(L, R);
        L = R + 1;
    }
    if (c_p1 > 0) {
        long long R = L + c_p1 - 1;
        lb += 1LL * sum_arith(L, R);
        L = R + 1;
    }
    // zeros contribute nothing; they occupy next c_0 positions but sum*0 == 0
    L += c_0;

    // negatives at the latest ranks, with -1 before -2, so that -2 get largest ranks
    long long tail_start = Ntot - (c_n1 + c_n2) + 1;
    if (c_n1 > 0) {
        long long R = tail_start + c_n1 - 1;
        lb += -1LL * sum_arith(tail_start, R);
    }
    if (c_n2 > 0) {
        long long R2_L = Ntot - c_n2 + 1;
        long long R2_R = Ntot;
        lb += -2LL * sum_arith(R2_L, R2_R);
    }
    // scale by dV outside (since g already includes dV).
    return lb * cfg_.dV;
}

long long IdaStarColumnPlacer::initial_bound(const State& s, int m, int h) const {
    (void)m; (void)h;
    // could also return heuristic_lb(s), but allowing bound to start at h(start)
    return s.g + heuristic_lb(s, m, h);
}

bool IdaStarColumnPlacer::ida_search(long long& bound, State& s,
                                     const int m, const int h,
                                     vector<pair<int,int>>& picks,
                                     IdaResult& result) {
    result.ida_iterations = 0;
    while (true) {
        ++result.ida_iterations;
        long long next_bound = std::numeric_limits<long long>::max();
        if (cfg_.verbose) {
            cout << "[IDA*] iter " << result.ida_iterations << " bound=" << bound << endl;
        }
        if (dfs(bound, next_bound, s, m, h, picks, result)) return true;
        if (next_bound == std::numeric_limits<long long>::max()) return false; // no path
        bound = next_bound;
    }
}

bool IdaStarColumnPlacer::dfs(long long bound, long long& next_bound, State& s,
                              const int m, const int h,
                              vector<pair<int,int>>& picks,
                              IdaResult& result) {
    ++result.dfs_calls;
    if (cfg_.max_dfs_calls && result.dfs_calls > cfg_.max_dfs_calls) return false;
    long long f = s.g + heuristic_lb(s, m, h);
    if (f > bound) {
        if (f < next_bound) next_bound = f;
        return false;
    }
        const int N = m*h;
        if (s.placed == N) {
        return true; // picks already hold the sequence
    }

    // generate legal moves (rows i) and order by Δcost ascending
    struct Move { int i, j; long long dcost; };
    vector<Move> moves;
    moves.reserve(m);
    for (int i = 0; i < m; ++i) {
        if (s.a[i] >= h) continue;
        if (i > 0 && s.a[i] + 1 > s.a[i-1]) continue; // maintain non-increasing frontier
        int j = s.a[i];
        long long rank = (long long)s.placed + 1;
        int w = weight_of(i, j, m, h);
        long long dcost = w * rank * cfg_.dV;
        moves.push_back({i, j, dcost});
    }
    if (moves.empty()) return false;

    std::sort(moves.begin(), moves.end(),
              [](const Move& a, const Move& b){ return a.dcost < b.dcost; });

    // expand
    for (const auto& mv : moves) {
        ++result.nodes_expanded;
        if (cfg_.max_nodes && result.nodes_expanded > cfg_.max_nodes) return false;

        // apply
        s.a[mv.i]++;
        s.placed++;
        s.g += mv.dcost;
        // update remaining weights multiset
        int w = weight_of(mv.i, mv.j, m, h);
        s.rem_cnt[w + 2]--; // index shift
        picks.emplace_back(mv.i, mv.j);

        // recurse
        if (dfs(bound, next_bound, s, m, h, picks, result)) {
            return true;
        }

        // undo
        picks.pop_back();
        s.rem_cnt[w + 2]++;
        s.g -= mv.dcost;
        s.placed--;
        s.a[mv.i]--;
    }
    return false;
}

IdaResult IdaStarColumnPlacer::solve(int m, int h) {
    IdaResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order.assign(m, std::vector<long long>(h, 0));
    R.x_of_col.assign(h, 0);

    // initial state
    State s;
    s.a.assign(m, 0);
    s.placed = 0;
    s.g = 0;

    // init remaining weight multiset
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int w = weight_of(i, j, m, h);
            s.rem_cnt[w + 2]++;
        }

    vector<pair<int,int>> picks;
    long long bound = initial_bound(s, m, h);
    bool ok = ida_search(bound, s, m, h, picks, R);

    if (!ok) {
        if (cfg_.verbose) {
            cout << "[IDA*] search failed (limits?)" << endl;
        }
        return R;
    }

    // reconstruct y_order from picks
    for (size_t k = 0; k < picks.size(); ++k) {
        int i = picks[k].first;
        int j = picks[k].second;
        R.y_order[i][j] = (long long)k + 1;
    }

    // verify
    R.oc_ok = check_global_OC(R.y_order);
    R.unique_ok = check_unique_1_to_N(R.y_order);
    R.objective_cost = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int w = weight_of(i, j, m, h);
            R.objective_cost += (long long)w * R.y_order[i][j] * cfg_.dV;
        }
    R.actual_hpwl = hpwl_edges_sum(R.y_order, R.x_of_col, cfg_.dV);

    if (cfg_.verbose) {
        cout << "[verify] objective=" << R.objective_cost
             << " actual HPWL=" << R.actual_hpwl
             << " OC=" << (R.oc_ok ? "OK" : "FAIL")
             << " Unique=" << (R.unique_ok ? "OK" : "FAIL") << endl;
        cout << "[stats] ida_iters=" << R.ida_iterations
             << " dfs_calls=" << R.dfs_calls
             << " nodes=" << R.nodes_expanded << endl;
    }
    return R;
}
