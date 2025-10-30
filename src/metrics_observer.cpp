#include "metrics_observer.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iomanip>

namespace metrics {

static inline void flatten(const Mat<long long>& A, std::vector<double>& out) {
    for (const auto& r : A) for (auto v : r) out.push_back((double)v);
}

static double median_of(std::vector<double> v) {
    if (v.empty()) return 0.0;
    size_t n = v.size();
    std::nth_element(v.begin(), v.begin() + n/2, v.end());
    double m = v[n/2];
    if ((n & 1) == 0) {
        auto it2 = std::max_element(v.begin(), v.begin()+n/2);
        m = (m + *it2) * 0.5;
    }
    return m;
}

static double gini_of_sorted(const std::vector<double>& x_sorted) {
    const size_t n = x_sorted.size();
    if (n == 0) return 0.0;
    double sum = std::accumulate(x_sorted.begin(), x_sorted.end(), 0.0);
    if (sum <= 0.0) return 0.0;
    // G = (2 * sum_{i=1..n} i*x_i)/(n*sum) - (n+1)/n
    double num = 0.0;
    for (size_t i = 0; i < n; ++i) num += (double)(i+1) * x_sorted[i];
    double G = (2.0 * num) / ( (double)n * sum ) - ( (double)n + 1.0 ) / (double)n;
    return G;
}

WeightStats compute_weight_stats(const Mat<long long>& WH, const Mat<long long>& WV, double top_ratio){
    std::vector<double> w; w.reserve((size_t)WH.size()*std::max(1,(int)WH[0].size()) + (size_t)WV.size()*std::max(1,(int)WV[0].size()));
    flatten(WH, w); flatten(WV, w);
    WeightStats S{};
    S.n = w.size();
    if (S.n == 0) return S;
    auto [minit, maxit] = std::minmax_element(w.begin(), w.end());
    S.minv = (long long)std::llround(*minit);
    S.maxv = (long long)std::llround(*maxit);
    double sum = std::accumulate(w.begin(), w.end(), 0.0);
    S.mean = sum / (double)S.n;
    double sq = 0.0;
    for (auto v : w) { double d=v-S.mean; sq += d*d; }
    S.stddev = std::sqrt(sq / (double)S.n);
    S.cv = (S.mean>0.0? S.stddev/S.mean : 0.0);

    // median
    std::vector<double> wcopy = w;
    S.median = median_of(wcopy);

    // gini & top share
    std::sort(w.begin(), w.end());
    S.gini = gini_of_sorted(w);
    size_t k = (size_t)std::ceil(top_ratio * (double)S.n);
    if (k < 1) k = 1;
    double top_sum = std::accumulate(w.end()-k, w.end(), 0.0);
    S.top_share = (sum>0.0? top_sum/sum : 0.0);
    return S;
}

EdgeStats compute_edge_stats(const Mat<int>& y, int m, int h, long long T){
    EdgeStats E{};
    // H edges
    for (int i=0;i<m;++i){
        for (int j=0;j<h-1;++j){
            long long L = std::llabs((long long)y[i][j+1] - (long long)y[i][j]);
            if (L > T) { E.viol_count++; E.viol_sum += (L - T); }
            if (L > E.maxL){
                E.maxL = L;
                E.maxEdge.horiz = true; E.maxEdge.i = i; E.maxEdge.j = j;
                E.maxEdge.cx = (double)i; E.maxEdge.cy = (double)j + 0.5;
                E.maxEdge.L = L; E.maxEdge.valid = true;
            }
        }
    }
    // V edges
    for (int i=0;i<m-1;++i){
        for (int j=0;j<h;++j){
            long long L = std::llabs((long long)y[i+1][j] - (long long)y[i][j]);
            if (L > T) { E.viol_count++; E.viol_sum += (L - T); }
            if (L > E.maxL){
                E.maxL = L;
                E.maxEdge.horiz = false; E.maxEdge.i = i; E.maxEdge.j = j;
                E.maxEdge.cx = (double)i + 0.5; E.maxEdge.cy = (double)j;
                E.maxEdge.L = L; E.maxEdge.valid = true;
            }
        }
    }
    return E;
}

JumpStats jump_between(const EdgeLoc& prev, const EdgeLoc& curr, bool has_prev){
    JumpStats J{};
    if (!has_prev || !prev.valid || !curr.valid) return J;
    double dx = curr.cx - prev.cx;
    double dy = curr.cy - prev.cy;
    J.manhattan = std::fabs(dx) + std::fabs(dy);
    J.euclid = std::sqrt(dx*dx + dy*dy);
    J.orient_flip = (prev.horiz != curr.horiz);
    return J;
}

void write_header(std::ostream& os){
    os << "# round T HPWL_eq HPWL_w WNS TNS maxL viol_cnt picked "
          "w_mean w_std w_cv w_median w_min w_max w_gini w_top10 "
          "jump_man jump_euc orient_flip maxEdge_type maxEdge_i maxEdge_j solve_ms feasible avg_jump\n";
}

void write_round(std::ostream& os, int round, long long T,
                 long long hpw_eq, long long hpw_w,
                 const EdgeStats& es, long long WNS, long long TNS, int picked,
                 const WeightStats& ws, const JumpStats& js,
                 double ms_ms, bool feasible, double avg_jump){
    os << round << " " << T << " " << hpw_eq << " " << hpw_w << " "
       << WNS << " " << TNS << " " << es.maxL << " " << es.viol_count << " " << picked << " "
       << std::fixed << std::setprecision(6)
       << ws.mean << " " << ws.stddev << " " << ws.cv << " " << ws.median << " "
       << (long long)ws.minv << " " << (long long)ws.maxv << " " << ws.gini << " " << ws.top_share << " "
       << js.manhattan << " " << js.euclid << " " << (js.orient_flip?1:0) << " "
       << (es.maxEdge.horiz? "H":"V") << " " << es.maxEdge.i << " " << es.maxEdge.j << " "
       << ms_ms << " " << (feasible?1:0) << " " << avg_jump << "\n";
}

} // namespace metrics
