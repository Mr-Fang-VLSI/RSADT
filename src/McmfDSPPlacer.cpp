#include "McmfDSPPlacer.h"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <unordered_set>

// LEMON
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>

using std::vector;
using std::pair;
using std::cout;
using std::endl;

static inline long long Tri(long long n) { return n*(n+1)/2; }

McmfDSPPlacer::McmfDSPPlacer(const Config& cfg) : cfg_(cfg) {}

// ===== R-SAD cost (LI) min over g =====
std::pair<long long,int> McmfDSPPlacer::LI_closed_form_int(int m, int w) const {
    if (w <= 0 || m <= 0) return {0, 1};
    if (w == 1) return { (long long)(m - 1), 0 }; // g ignored
    int g_ub = std::min(m, w) / 2;
    long long best6 = std::numeric_limits<long long>::max();
    int bestg = 1;
    for (int g = 1; g <= g_ub; ++g) {
        long long h = w;
        long long g2 = 1LL*g*g, g3 = g2*g;
        long long val6 = -4*g3 + 12*h*g2 + (4 - 6*h*h - 6*h)*1LL*g
                       + 6LL*m*h*h + 6LL*m*h - 6LL*m - 6LL*h;
        if (val6 < best6) { best6 = val6; bestg = g; }
    }
    return {best6/6, bestg};
}

// ===== LI(m,w,g) (no minimization) =====
long long McmfDSPPlacer::LI_by_g_int(int m, int w, int g) const {
    if (w <= 0 || m <= 0) return 0;
    if (w == 1) return (long long)(m - 1); // g ignored
    long long h = w;
    long long g2 = 1LL*g*g, g3 = g2*g;
    long long val6 = -4*g3 + 12*h*g2 + (4 - 6*h*h - 6*h)*1LL*g
                   + 6LL*m*h*h + 6LL*m*h - 6LL*m - 6LL*h;
    return val6/6;
}

// ===== RSAD region index helpers =====
long long McmfDSPPlacer::OLL(int i, int j) {
    if (i >= j) return 1LL*i*i - i + j;
    return 1LL*(j-1)*(j-1) + i;
}
long long McmfDSPPlacer::OLR(int i, int j, int g) {
    if (i <= g - j + 1) {
        long long num = -1LL*j*j + (2LL*g + 3)*j - 2LL*g + 2LL*i - 2;
        return num / 2;
    } else {
        return Tri(g-1) + Tri(i-1) + j;
    }
}
long long McmfDSPPlacer::OUL(int i, int j, int g) {
    if (i <= g - j + 1) {
        long long num = -1LL*i*i + (2LL*g + 3)*i - 2LL*g + 2LL*j - 2;
        return num / 2;
    } else {
        return Tri(g-1) + Tri(j-1) + i;
    }
}
long long McmfDSPPlacer::OUR(int i, int j, int g) {
    return 1LL*g*g + 1 - OLL(g - i + 1, g - j + 1);
}

// ===== RSAD mapping (7 regions) =====
std::vector<std::vector<long long>>
McmfDSPPlacer::rsad_mapping(int m, int w, int g) const {
    vector<vector<long long>> y(m, vector<long long>(w, 0));
    if (w == 1) { for (int i = 0; i < m; ++i) y[i][0] = i + 1; return y; }

    int mid_cols = std::max(0, w - 2*g);
    int mid_rows = std::max(0, m - 2*g);
    long long offset = 0;

    // 1) Lower-left g x g
    for (int i = 1; i <= g; ++i)
    for (int j = 1; j <= g; ++j)
        y[i-1][j-1] = offset + OLL(i, j);
    offset += 1LL*g*g;

    // 2) Lower-middle: g x (w-2g)
    for (int j = g+1; j <= w - g; ++j) {
        int jrel = j - g;
        for (int i = 1; i <= g; ++i)
            y[i-1][j-1] = offset + 1LL*(jrel - 1)*g + i;
    }
    offset += 1LL*mid_cols*g;

    // 3) Lower-right g x g
    for (int j = w - g + 1; j <= w; ++j) {
        int jrel = j - (w - g);
        for (int i = 1; i <= g; ++i)
            y[i-1][j-1] = offset + OLR(i, jrel, g);
    }
    offset += 1LL*g*g;

    // 4) Central band: (m-2g) x w
    for (int i = g + 1; i <= m - g; ++i) {
        int irel = i - g;
        for (int j = 1; j <= w; ++j)
            y[i-1][j-1] = offset + 1LL*(irel - 1)*w + j;
    }
    offset += 1LL*mid_rows*w;

    // 5) Upper-left g x g
    for (int i = m - g + 1; i <= m; ++i) {
        int irel = i - (m - g);
        for (int j = 1; j <= g; ++j)
            y[i-1][j-1] = offset + OUL(irel, j, g);
    }
    offset += 1LL*g*g;

    // 6) Upper-middle: g x (w-2g)
    for (int j = g+1; j <= w - g; ++j) {
        int jrel = j - g;
        for (int i = m - g + 1; i <= m; ++i) {
            int irel = i - (m - g);
            y[i-1][j-1] = offset + 1LL*(jrel - 1)*g + irel;
        }
    }
    offset += 1LL*mid_cols*g;

    // 7) Upper-right g x g
    for (int j = w - g + 1; j <= w; ++j) {
        int jrel = j - (w - g);
        for (int i = m - g + 1; i <= m; ++i) {
            int irel = i - (m - g);
            y[i-1][j-1] = offset + OUR(irel, jrel, g);
        }
    }
    offset += 1LL*g*g;

    assert(offset == 1LL*m*w && "RSAD numbering must cover 1..m*w");
    return y;
}

// ===== Vpen between patterns (for cost) =====
long long McmfDSPPlacer::vpen_between_patterns(int m, int w, int g_prev, int g_next) const {
    auto y_prev = rsad_mapping(m, w, (w==1?1:g_prev));
    auto y_next = rsad_mapping(m, w, (w==1?1:g_next));
    long long s_needed = 0;
    long long sumL = 0, sumR = 0;
    for (int i = 0; i < m; ++i) {
        long long R = y_prev[i][w-1];
        long long L = y_next[i][0];
        s_needed = std::max(s_needed, R - L);
        sumR += R; sumL += L;
    }
    return 1LL*m*s_needed + (sumL - sumR); // steps; multiply dV outside
}

// ===== max-T helpers =====
long long McmfDSPPlacer::max_intra_step(int m, int w, int g) const {
    if (w <= 1) return 0;
    auto y = rsad_mapping(m, w, (w==1?1:g));
    long long mx = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j+1 < w; ++j)
            mx = std::max(mx, y[i][j+1] - y[i][j]); // nonnegative under OC
    return mx;
}
long long McmfDSPPlacer::max_inter_step_after_shift(int m, int w, int g_prev, int g_next) const {
    auto yp = rsad_mapping(m, w, (w==1?1:g_prev));
    auto yn = rsad_mapping(m, w, (w==1?1:g_next));
    long long max_rl = std::numeric_limits<long long>::min();
    long long min_rl = std::numeric_limits<long long>::max();
    for (int i = 0; i < m; ++i) {
        long long diff = yp[i][w-1] - yn[i][0]; // R-L
        max_rl = std::max(max_rl, diff);
        min_rl = std::min(min_rl, diff);
    }
    return max_rl - min_rl; // = max after minimal shift
}

// ===== HPWL (grid edges) =====
long long McmfDSPPlacer::hpwl_edges_sum(const vector<vector<long long>>& y,
                                        const vector<long long>& x,
                                        long long dV) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long S = 0;
    // horizontal neighbors
    for (int i = 0; i < m; ++i)
        for (int j = 0; j+1 < h; ++j) {
            long long dx = std::llabs(x[j+1] - x[j]);
            long long dy = std::llabs(y[i][j+1] - y[i][j]) * dV;
            S += dx + dy;
        }
    // vertical neighbors
    for (int i = 0; i+1 < m; ++i)
        for (int j = 0; j < h; ++j) {
            long long dy = std::llabs(y[i+1][j] - y[i][j]) * dV;
            S += dy;
        }
    return S;
}

// ===== global OC check =====
bool McmfDSPPlacer::check_global_OC(const vector<vector<long long>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
    for (int j1 = 0; j1 < h; ++j1)
    for (int i2 = i1; i2 < m; ++i2)
    for (int j2 = j1; j2 < h; ++j2)
        if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

long long McmfDSPPlacer::max_horizontal_step_overall(const vector<vector<long long>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long mx = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j+1 < h; ++j)
            mx = std::max(mx, std::llabs(y[i][j+1] - y[i][j]));
    return mx;
}

// ===== MCMF: equal width, OC-aware, with max-T filters =====
std::tuple<int,int,long long,std::vector<int>>
McmfDSPPlacer::solve_mcmf_equal_width_oc(int m, int h) const {
    using Graph = lemon::ListDigraph;

    // candidate widths: divisors of h, with s=h/w <= max_site_cols
    vector<int> W;
    for (int w = 1; w <= h; ++w)
        if (h % w == 0 && (h / w) <= cfg_.max_site_cols) W.push_back(w);
    if (W.empty()) W.push_back(h);

    long long best_cost = std::numeric_limits<long long>::max();
    int best_w = W.front(), best_s = h / best_w;
    std::vector<int> best_g_list;

    for (int w : W) {
        int s = h / w;

        // allowed g for this w after intra-step check
        vector<int> G_all;
        if (w == 1) G_all.push_back(0);
        else {
            int g_ub = std::min(m, w) / 2;
            for (int g = 1; g <= g_ub; ++g) {
                if (cfg_.max_T_steps > 0) {
                    long long mi = max_intra_step(m, w, g);
                    if (mi > cfg_.max_T_steps) continue; // violates max-T inside segment
                }
                G_all.push_back(g);
            }
        }
        if (G_all.empty()) continue;

        Graph g;
        auto cap    = Graph::ArcMap<int>(g);
        auto cost   = Graph::ArcMap<long long>(g);
        auto supply = Graph::NodeMap<long long>(g);

        Graph::Node src = g.addNode();
        Graph::Node sink = g.addNode();
        supply[src] = +1; supply[sink] = -1;

        // layer nodes N[t][ki] , t=0..s-1
        std::vector<std::vector<Graph::Node>> N(s);
        for (int t = 0; t < s; ++t) {
            N[t].resize(G_all.size());
            for (size_t ki = 0; ki < G_all.size(); ++ki) {
                N[t][ki] = g.addNode();
                supply[N[t][ki]] = 0;
            }
        }

        // src arcs (cost = LI of first seg)
        std::vector<Graph::Arc> src_arcs(G_all.size(), lemon::INVALID);
        for (size_t ki = 0; ki < G_all.size(); ++ki) {
            int g1 = G_all[ki];
            auto a = g.addArc(src, N[0][ki]);
            cap[a] = 1;
            long long c = LI_by_g_int(m, w, g1) * cfg_.dV;
            cost[a] = c;
            src_arcs[ki] = a;
        }

        // transitions with max-T filtering on boundary
        std::vector<std::vector<std::vector<Graph::Arc>>> trans_arcs; // [t][ki][kj]
        trans_arcs.resize(std::max(0, s-1));
        for (int t = 0; t+1 < s; ++t) {
            trans_arcs[t].resize(G_all.size());
            for (size_t ki = 0; ki < G_all.size(); ++ki) {
                trans_arcs[t][ki].resize(G_all.size(), lemon::INVALID);
                int g_prev = G_all[ki];
                for (size_t kj = 0; kj < G_all.size(); ++kj) {
                    int g_next = G_all[kj];
                    if (cfg_.max_T_steps > 0) {
                        long long mx = max_inter_step_after_shift(m, w, g_prev, g_next);
                        if (mx > cfg_.max_T_steps) continue; // forbid this transition
                    }
                    auto a = g.addArc(N[t][ki], N[t+1][kj]);
                    cap[a] = 1;
                    long long vpen = cfg_.enforce_global_oc ? vpen_between_patterns(m, w, g_prev, g_next) : 0;
                    long long c = 1LL*m*cfg_.dH + vpen*cfg_.dV + LI_by_g_int(m, w, g_next)*cfg_.dV;
                    cost[a] = c;
                    trans_arcs[t][ki][kj] = a;
                }
            }
        }

        // last layer -> sink
        std::vector<Graph::Arc> sink_arcs(G_all.size(), lemon::INVALID);
        for (size_t ki = 0; ki < G_all.size(); ++ki) {
            auto a = g.addArc(N[s-1][ki], sink);
            cap[a] = 1; cost[a] = 0;
            sink_arcs[ki] = a;
        }

        lemon::NetworkSimplex<Graph,int,long long> ns(g);
        ns.upperMap(cap).costMap(cost).supplyMap(supply);
        auto status = ns.run();
        if (status != decltype(ns)::OPTIMAL) continue; // no feasible path under max-T

        long long Z = ns.totalCost();

        // extract path
        int pick0 = -1;
        for (size_t ki = 0; ki < G_all.size(); ++ki)
            if (src_arcs[ki] != lemon::INVALID && ns.flow(src_arcs[ki]) > 0) { pick0 = (int)ki; break; }
        if (pick0 < 0) continue;

        std::vector<int> g_list; g_list.reserve(s);
        g_list.push_back(G_all[pick0]);
        int cur_idx = pick0;
        for (int t = 0; t+1 < s; ++t) {
            int nxt = -1;
            for (size_t kj = 0; kj < G_all.size(); ++kj) {
                auto a = trans_arcs[t][cur_idx][kj];
                if (a != lemon::INVALID && ns.flow(a) > 0) { nxt = (int)kj; break; }
            }
            if (nxt < 0) break;
            g_list.push_back(G_all[nxt]);
            cur_idx = nxt;
        }
        if ((int)g_list.size() != s) continue;

        if (Z < best_cost) {
            best_cost = Z; best_w = w; best_s = s; best_g_list = std::move(g_list);
        }
    }

    return {best_w, best_s, best_cost, best_g_list};
}

// ===== Public solve =====
PlacementResult McmfDSPPlacer::solve(int m, int h) {
    PlacementResult R;
    R.m = m; R.h = h; R.dH = cfg_.dH; R.dV = cfg_.dV;
    R.max_T_limit = cfg_.max_T_steps;

    // 1) OC-aware MCMF with equal width + max-T
    int w = 0, s = 0; long long Z = 0; std::vector<int> glist;
    std::tie(w, s, Z, glist) = solve_mcmf_equal_width_oc(m, h);
    R.mcmf_cost = Z;

    if (cfg_.verbose) {
        cout << "[MCMF] h=" << h << " OC-aware equal width"
             << (cfg_.max_T_steps>0? (" with maxT="+std::to_string(cfg_.max_T_steps)) : "")
             << " -> w=" << w << ", s=" << s << ", cost=" << Z << endl;
    }

    // 2) materialize mapping with baseline shifts to satisfy global OC
    R.segments.clear();
    R.segments.reserve(s);
    R.x_of_col.assign(h, 0);
    R.y_order.assign(m, vector<long long>(h, 0));
    R.y_label.assign(m, vector<long long>(h, 0));

    long long base = 0;
    long long intra_sum = 0, interH_sum = 1LL*(s-1)*m*cfg_.dH, interV_sum = 0;
    vector<long long> prevR(m, 0);
    int col0 = 0;

    for (int t = 0; t < s; ++t) {
        int gstar = (w==1 ? 0 : glist[t]);
        auto yseg = rsad_mapping(m, w, (w==1?1:gstar));

        if (t > 0) {
            long long sneed = 0, sumL = 0, sumR = 0;
            for (int i = 0; i < m; ++i) {
                long long L = yseg[i][0];
                sneed = std::max(sneed, prevR[i] - L);
                sumL += L; sumR += prevR[i];
            }
            base += sneed;
            interV_sum += (1LL*m*sneed + (sumL - sumR)) * cfg_.dV;
        }

        long long label_base = 1LL*t * m * w; // global 1..(m*h)
        for (int j = 0; j < w; ++j) {
            R.x_of_col[col0 + j] = 1LL*t * cfg_.dH;
            for (int i = 0; i < m; ++i) {
                R.y_order[i][col0 + j] = yseg[i][j] + base;
                R.y_label[i][col0 + j] = yseg[i][j] + label_base;
            }
        }
        for (int i = 0; i < m; ++i) prevR[i] = yseg[i][w-1] + base;

        Segment seg; seg.start_col = col0; seg.width = w; seg.g = (w==1?0:gstar); seg.mirrored = false;
        R.segments.push_back(seg);

        col0 += w;
        intra_sum += LI_by_g_int(m, w, (w==1?0:gstar)) * cfg_.dV;
    }
    R.used_site_cols = s;

    // 3) Actual HPWL
    R.actual_hpwl = hpwl_edges_sum(R.y_order, R.x_of_col, cfg_.dV);

    // 4) Checks
    R.oc_ok = check_global_OC(R.y_order);
    {
        std::unordered_set<unsigned long long> occ;
        bool ok = true;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < h; ++j) {
                unsigned long long x = (unsigned long long)(R.x_of_col[j] & 0xffffffffULL);
                unsigned long long y = (unsigned long long)(R.y_order[i][j] & 0xffffffffULL);
                unsigned long long key = (x << 32) ^ y;
                if (!occ.insert(key).second) ok = false;
            }
        R.unique_ok = ok;
    }
    R.max_step_observed = max_horizontal_step_overall(R.y_order);
    R.maxT_ok = (cfg_.max_T_steps <= 0) ? true : (R.max_step_observed <= cfg_.max_T_steps);

    if (cfg_.verbose) {
        cout << "  [summary] intra=" << intra_sum
             << ", interH=" << interH_sum
             << ", interV=" << interV_sum
             << ", total=" << (intra_sum + interH_sum + interV_sum)
             << " (should equal MCMF)" << endl;
        cout << "  [verify] MCMF cost=" << R.mcmf_cost
             << ", actual HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok ? "OK" : "FAIL")
             << ", Unique=" << (R.unique_ok ? "OK" : "FAIL");
        if (cfg_.max_T_steps > 0) {
            cout << ", maxT(limit=" << cfg_.max_T_steps << ", observed=" << R.max_step_observed
                 << ")=" << (R.maxT_ok ? "OK" : "FAIL");
        }
        cout << endl;
    }
    return R;
}
