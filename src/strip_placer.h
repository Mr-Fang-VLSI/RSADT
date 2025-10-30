#pragma once
#include "lightOCShortest.h"
#include "momentum_weighter.h"
#include <vector>
#include <string>
#include <limits>
#include <iostream>

template<typename T> using Mat = std::vector<std::vector<T>>;

class StripPlacer {
public:
    struct Params {
        int m = 0, h = 0, K = 2;          // 切成 K 份（仅支持 2 或 4），h 必须能被 K 整除
        int rounds = 5;                   // 外层权重迭代轮数
        long long dV = 1;                 // φ缩放
        long long T = 0;                  // Max-T
        double alpha = 0.5, eta = 0.5, top = 0.10, cap = 32.0, scale = 1000.0;
        bool progress = true;             // 打印
        bool use_pow2 = true;             // 2^k bitpack

        // 可选：直接传入条带尺寸 m x (w-1)/(m-1 x w) 的权重（如不传则内部全1）
        Mat<long long> wH_local, wV_local;
        bool has_local_weights = false;
    };

    explicit StripPlacer(const Params& P) : P_(P) {
        if(P_.K!=2 && P_.K!=4) throw std::runtime_error("[StripPlacer] K must be 2 or 4.");
        if(P_.h % P_.K != 0)   throw std::runtime_error("[StripPlacer] h must be divisible by K.");
        w_ = P_.h / P_.K;
        if(P_.T<=0) P_.T = std::max(P_.m, w_); // 合理默认
    }

    // 在单个条带上求解一次（m x w），把最优 y_local_ 保存（1-based）
    void solve_once() {
        const int m=P_.m, w=w_;
        // 1) 准备权重（整条带等权=scale；若给了 wH/wV_local 就用之）
        Mat<long long> wH(m, std::vector<long long>(w-1, (long long)std::llround(P_.scale)));
        Mat<long long> wV(m-1, std::vector<long long>(w,   (long long)std::llround(P_.scale)));
        if(P_.has_local_weights){
            if((int)P_.wH_local.size()==m && (w>0 ? (int)P_.wH_local[0].size()==w-1 : true)) wH = P_.wH_local;
            if((int)P_.wV_local.size()==m-1 && (w>0 ? (int)P_.wV_local[0].size()==w   : true)) wV = P_.wV_local;
        }
        MomentumWeighter::Params WP; WP.alpha=P_.alpha; WP.eta=P_.eta; WP.top_ratio=P_.top; WP.cap_max=P_.cap; WP.export_scale=P_.scale;
        MomentumWeighter weighter(m,w,WP);
        weighter.set_initial(wH, wV, P_.scale);

        long long best_hp = std::numeric_limits<long long>::max();
        Mat<int>    best_y;

        for(int it=1; it<=P_.rounds; ++it){
            // 导出整数权重
            Mat<long long> WH_ll, WV_ll; weighter.export_ll(WH_ll, WV_ll);

            lightOCShortest::Config cfg;
            cfg.dV = (P_.dV==0?1:P_.dV);  // 防御：避免 dV=0 退化
            cfg.verbose  = P_.progress;
            cfg.progress = P_.progress;
            cfg.max_steps   = -1;
            cfg.omp_threads = 0;
            cfg.low_mem     = true;
            cfg.use_pow2_pack = P_.use_pow2;
            cfg.W.enabled=true; cfg.W.wH=WH_ll; cfg.W.wV=WV_ll;

            lightOCShortest solver(cfg);
            auto R = solver.solve(m, w);

            long long hpW = lightOCShortest::hpwl_sum_weighted(R.y_order, cfg.W, 1);
            long long hpE = lightOCShortest::hpwl_sum_equal   (R.y_order, 1);
            if(P_.progress){
                std::cout<<"  [Strip] Round "<<it<<"  HPWL_w="<<hpW<<"  HPWL_eq="<<hpE<<"  report="<<R.total_cost<<"\n";
            }
            if(hpW < best_hp){ best_hp=hpW; best_y=R.y_order; }

            auto rep = weighter.update_from_layout(R.y_order, P_.T, /*pick_by_ratio=*/true);
            if(P_.progress){
                std::cout<<"    [Updt] WNS="<<rep.WNS<<" TNS="<<rep.TNS<<" maxL="<<rep.max_edge_len<<" picked="<<rep.picked_edges<<"\n";
            }
            if(rep.WNS >= 0) break; // 满足 T，提前收敛
        }

        if(best_y.empty()) throw std::runtime_error("[StripPlacer] empty solution.");
        y_local_ = std::move(best_y); // 1-based
    }

    // 构建全局 X/Y（m x h）；复制列内顺序到 K 个列，按 j->(col=j/w_, jloc=j%w_)
    void build_global_xy(Mat<int>& X, Mat<int>& Y) const {
        if(y_local_.empty()) throw std::runtime_error("[StripPlacer] call solve_once() first.");
        const int m=P_.m, h=P_.h, w=w_;
        X.assign(m, std::vector<int>(h, 0));
        Y.assign(m, std::vector<int>(h, 0));
        for(int i=0;i<m;++i){
            for(int j=0;j<h;++j){
                int col  = j / w;           // 0..K-1
                int jloc = j % w;           // 0..w-1
                int y1b  = y_local_[i][jloc];   // 1-based
                X[i][j] = col;
                Y[i][j] = y1b - 1;          // 0-based 输出
            }
        }
    }

    int w() const { return w_; }

private:
    Params P_;
    int w_ = 0;
    Mat<int> y_local_; // m x w，1-based
};
