// src/momentum_weighter.h
#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <limits>
#include <tuple>
#include <string>
#include <cassert>

struct MomentumWeighterReport {
    // ????(????/??)
    long long HPWL_equal = 0;   // ?? HPWL(????)
    long long HPWL_weighted = 0;
    long long max_edge_len = 0; // ??????
    long long WNS = 0;          // ?? slack = min(T - L_e)
    long long TNS = 0;          // ? slack ??(?????)
    int       picked_edges = 0; // ???????????
};

// ????
template<typename T>
using Mat = std::vector<std::vector<T>>;

// —— ????? ——
//
// ??:
// m×h ??;H ??? = m × (h-1);V ??? = (m-1) × h?
// “??”??? double ?,???????? long long (??? scale)?
//
// ????(???,?? DREAMPlace 4.0 ?(9)-(12)):
//   dlogw^(m+1) = a * dlogw^(m) + ? * log(1 + c_mom^(m))
//   logw^(m+1)  = logw^(m)     + dlogw^(m+1)
// ?? c_mom ?? Max-T ??:
//   s_e = T - L_e,WNS = min_e s_e
//   ? WNS >= 0,? c_mom = 0;?? c_mom = max( s_e / WNS, 0 )?
// ?? Top-p% ?????(??),????????? a?
//
class MomentumWeighter {
public:
    struct Params {
        double alpha    = 0.5;   // ???? a ? [0,1]
        double eta      = 0.5;   // ?? ? > 0;???? ?=1-a
        double top_ratio= 0.10;  // ????? Top ??(? slack ??)
        double cap_max  = 32.0;  // ????(??);<=0 ?????
        double cap_min  = 1.0;   // ????
        double export_scale = 1000.0; // ?????????
    };

    MomentumWeighter(int m, int h, const Params& P)
        : m_(m), h_(h), P_(P)
        , logwH_(m, std::vector<double>(h-1, 0.0))
        , logwV_(m-1, std::vector<double>(h, 0.0))
        , dlogH_(m, std::vector<double>(h-1, 0.0))
        , dlogV_(m-1, std::vector<double>(h, 0.0))
    {
        // ????=1 ? log(1)=0
    }

    // ????????????(??????? wH,wV;?????????)
    void set_initial(const Mat<long long>& wH_ll, const Mat<long long>& wV_ll, double scale){
        assert((int)wH_ll.size()==m_ && (int)wH_ll[0].size()==h_-1);
        assert((int)wV_ll.size()==m_-1 && (int)wV_ll[0].size()==h_);
        for(int i=0;i<m_;++i) for(int j=0;j<h_-1;++j){
            double w = std::max(1.0, (double)wH_ll[i][j]/scale);
            logwH_[i][j] = std::log(w);
        }
        for(int i=0;i<m_-1;++i) for(int j=0;j<h_;++j){
            double w = std::max(1.0, (double)wV_ll[i][j]/scale);
            logwV_[i][j] = std::log(w);
        }
    }

    // ??? double ????? long long(??????)
    void export_ll(Mat<long long>& wH_ll, Mat<long long>& wV_ll) const {
        wH_ll.assign(m_, std::vector<long long>(h_-1, 0));
        wV_ll.assign(m_-1, std::vector<long long>(h_, 0));
        for(int i=0;i<m_;++i) for(int j=0;j<h_-1;++j){
            double w = std::exp(logwH_[i][j]);
            if(P_.cap_max>0) w = std::min(w, P_.cap_max);
            w = std::max(w, P_.cap_min);
            long long x = (long long)std::llround(w * P_.export_scale);
            if(x<=0) x=1;
            wH_ll[i][j]=x;
        }
        for(int i=0;i<m_-1;++i) for(int j=0;j<h_;++j){
            double w = std::exp(logwV_[i][j]);
            if(P_.cap_max>0) w = std::min(w, P_.cap_max);
            w = std::max(w, P_.cap_min);
            long long x = (long long)std::llround(w * P_.export_scale);
            if(x<=0) x=1;
            wV_ll[i][j]=x;
        }
    }

    // ??:???? layout(y) ? T ???????,??????
    MomentumWeighterReport update_from_layout(
        const Mat<int>& y, long long T,
        bool pick_by_ratio = true /*true:? c_mom ?;false:? slack ?*/
    ){
        MomentumWeighterReport R{};
        const int m=m_, h=h_;

        // ???? & slack
        struct EdgeRef { bool horiz; int i,j; long long L; long long s; double ratio; };
        std::vector<EdgeRef> bad; bad.reserve((size_t)m*(h-1) + (size_t)(m-1)*h);

        R.max_edge_len = 0;
        auto add_edge = [&](bool horiz, int i,int j, long long L){
            long long s = T - L;            // slack
            R.max_edge_len = std::max(R.max_edge_len, L);
            if(s<0) R.TNS += -s;
            // ratio = max(s/WNS, 0) ?????,?? WNS ???
            bad.push_back({horiz,i,j,L,s,0.0});
        };

        // ???? HPWL (???/??) & ??? slack
        R.HPWL_equal = 0;
        for(int i=0;i<m;++i){
            for(int j=0;j<h-1;++j){
                long long L = std::llabs((long long)y[i][j+1] - (long long)y[i][j]);
                R.HPWL_equal += L;
                add_edge(true, i,j,L);
            }
        }
        for(int i=0;i<m-1;++i){
            for(int j=0;j<h;++j){
                long long L = std::llabs((long long)y[i+1][j] - (long long)y[i][j]);
                R.HPWL_equal += L;
                add_edge(false, i,j,L);
            }
        }

        // WNS
        R.WNS = std::numeric_limits<long long>::max();
        for(auto& e: bad) R.WNS = std::min(R.WNS, e.s);

        // ?? ratio ???“??”
        std::vector<EdgeRef> cand; cand.reserve(bad.size());
        if(R.WNS < 0){
            for(auto &e: bad){
                double r = (double)e.s / (double)R.WNS; // WNS<0 ? ?/?=?
                if(r > 0.0) { e.ratio = r; cand.push_back(e); }
            }
        }
        // Top ??
        int K = (int)std::llround((double)cand.size() * P_.top_ratio);
        if(K < 1 && !cand.empty()) K = 1;
        if(K > (int)cand.size())   K = (int)cand.size();

        if(K>0){
            if(pick_by_ratio){
                std::nth_element(cand.begin(), cand.begin()+K-1, cand.end(),
                                 [](const EdgeRef& a, const EdgeRef& b){ return a.ratio > b.ratio; });
            }else{
                std::nth_element(cand.begin(), cand.begin()+K-1, cand.end(),
                                 [](const EdgeRef& a, const EdgeRef& b){ return a.s < b.s; }); // ????
            }
            cand.resize(K);
        }else{
            cand.clear();
        }
        R.picked_edges = (int)cand.size();

        // —— ???? —— //
        // ??????????(??? Top ???? dlog *= a)
        for(int i=0;i<m_;++i) for(int j=0;j<h_-1;++j) dlogH_[i][j] *= P_.alpha;
        for(int i=0;i<m_-1;++i) for(int j=0;j<h_;  ++j) dlogV_[i][j] *= P_.alpha;

        // ?? Top ??: dlog = a*dlog + ?*log(1+c_mom)
        auto bump = [&](bool horiz, int i,int j, double cmom){
            double inc = P_.eta * std::log1p(cmom);      // log(1 + c)
            if(horiz) dlogH_[i][j] += inc;
            else      dlogV_[i][j] += inc;
        };
        if(R.WNS < 0){
            for(const auto& e: cand){
                double cmom = std::max(0.0, (double)e.s / (double)R.WNS);
                if(cmom>0.0) bump(e.horiz, e.i, e.j, cmom);
            }
        }
        // ??? logw
        for(int i=0;i<m_;++i) for(int j=0;j<h_-1;++j) logwH_[i][j] += dlogH_[i][j];
        for(int i=0;i<m_-1;++i) for(int j=0;j<h_;  ++j) logwV_[i][j] += dlogV_[i][j];

        return R;
    }

private:
    int m_, h_;
    Params P_;

    // log ?? & ??
    Mat<double> logwH_, logwV_;
    Mat<double> dlogH_, dlogV_;
};
