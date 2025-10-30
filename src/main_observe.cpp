#include "lightOCShortest.h"
#include "momentum_weighter.h"
#include "metrics_observer.h"

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <cctype>
#include <limits>

using std::cout; using std::cerr; using std::endl; using std::string;
template<typename T> using Mat = std::vector<std::vector<T>>;

static bool parse_flag(const string& s, const string& key, string& out){
    auto p=s.find('='); if(p==string::npos) return false;
    if(s.substr(0,p)==key){ out=s.substr(p+1); return true; }
    return false;
}
static bool is_integer(const char* s){
    if(!s||!*s) return false; const char* p=s; if(*p=='+'||*p=='-') ++p; if(!*p) return false;
    while(*p){ if(!std::isdigit((unsigned char)*p)) return false; ++p; } return true;
}

int main(int argc,char** argv){
    if(argc<3){
        cout<<"Usage: ./build/oc_observe m h [--rounds=30] [--T=auto] [--dV=1]"
              " [--alpha=0.5] [--eta=0.5] [--top=0.10] [--cap=32] [--scale=1000]"
              " [--pow2=1] [--progress=1] [--outfile path]\n";
        return 0;
    }
    const int m = std::atoi(argv[1]);
    const int h = std::atoi(argv[2]);

    // defaults
    int rounds = 30;
    long long T = (long long)std::llround(1.5 * (double)std::min(m,h)); // 典型 T
    long long dV = 1;
    double alpha=0.5, eta=0.5, top=0.10, cap=32.0, scale=1000.0;
    bool use_pow2=true, progress=true;
    string out_path;

    string v;
    for(int i=1;i<argc;++i){
        string s(argv[i]);
        if(parse_flag(s,"--rounds",v)) rounds=std::stoi(v);
        else if(parse_flag(s,"--T",v)) T=std::llround(std::stod(v));
        else if(parse_flag(s,"--dV",v)) dV=std::llround(std::stod(v));
        else if(parse_flag(s,"--alpha",v)) alpha=std::stod(v);
        else if(parse_flag(s,"--eta",v))   eta=std::stod(v);
        else if(parse_flag(s,"--top",v))   top=std::stod(v);
        else if(parse_flag(s,"--cap",v))   cap=std::stod(v);
        else if(parse_flag(s,"--scale",v)) scale=std::stod(v);
        else if(parse_flag(s,"--pow2",v))  use_pow2=(v!="0");
        else if(parse_flag(s,"--progress",v)) progress=(v=="1");
        else if(s=="--outfile" && i+1<argc) out_path=argv[++i];
    }
    if(out_path.empty()){
        out_path = "metrics_m"+std::to_string(m)+"_h"+std::to_string(h)+"_T"+std::to_string(T)+".txt";
    }

    // initialize unit weights (scale) and weighter
    Mat<long long> wH(m, std::vector<long long>(h-1, (long long)std::llround(scale)));
    Mat<long long> wV(m-1, std::vector<long long>(h,   (long long)std::llround(scale)));
    MomentumWeighter::Params WP; WP.alpha=alpha; WP.eta=eta; WP.top_ratio=top; WP.cap_max=cap; WP.export_scale=scale;
    MomentumWeighter weighter(m,h,WP);
    weighter.set_initial(wH, wV, scale);

    cout<<"[Observe] m="<<m<<" h="<<h<<" rounds="<<rounds<<" T="<<T
        <<" alpha="<<alpha<<" eta="<<eta<<" top="<<top<<" cap="<<cap
        <<" pow2="<<(use_pow2?1:0)<<endl;

    std::ofstream fout(out_path);
    if(!fout){ cerr<<"[Observe] cannot open "<<out_path<<endl; return 2; }
    metrics::write_header(fout);

    // prev max-edge for jump metric
    EdgeLoc prev_edge; bool has_prev=false;
    double sum_jump=0.0; int jump_cnt=0;

    Mat<int> y_best; long long best_hp=std::numeric_limits<long long>::max();

    for(int it=1; it<=rounds; ++it){
        // export integer weights to DP
        Mat<long long> WH_ll, WV_ll; weighter.export_ll(WH_ll, WV_ll);

        // weight stats BEFORE solve (本轮进入 DP 的权重)
        auto wstats = metrics::compute_weight_stats(WH_ll, WV_ll, 0.10);

        // config
        lightOCShortest::Config cfg;
        cfg.dV=(dV==0?1:dV); cfg.verbose=progress; cfg.progress=progress; cfg.max_steps=-1;
        cfg.omp_threads=0; cfg.low_mem=true; cfg.use_pow2_pack=use_pow2;
        cfg.W.enabled=true; cfg.W.wH=WH_ll; cfg.W.wV=WV_ll;

        // solve
        lightOCShortest solver(cfg);
        auto t0=std::chrono::high_resolution_clock::now();
        auto R = solver.solve(m,h);
        auto t1=std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double,std::milli>(t1-t0).count();

        long long hpW = lightOCShortest::hpwl_sum_weighted(R.y_order, cfg.W, 1);
        long long hpE = lightOCShortest::hpwl_sum_equal   (R.y_order, 1);

        // report from weighter (WNS/TNS/maxL/picked)
        auto rep = weighter.update_from_layout(R.y_order, T, /*pick_by_ratio=*/true); // 定义见源码注释 :contentReference[oaicite:2]{index=2}

        // our edge stats for max-edge location & viol_count
        auto estat = metrics::compute_edge_stats(R.y_order, m, h, T);

        // jump metrics
        auto jstat = metrics::jump_between(prev_edge, estat.maxEdge, has_prev);
        if(has_prev){ sum_jump += jstat.manhattan; ++jump_cnt; }
        double avg_jump = (jump_cnt>0? sum_jump/(double)jump_cnt : 0.0);
        prev_edge = estat.maxEdge; has_prev=true;

        bool feasible = (rep.WNS >= 0);

        // record
        metrics::write_round(fout, it, T, hpE, hpW, estat, rep.WNS, rep.TNS, rep.picked_edges,
                             wstats, jstat, ms, feasible, avg_jump);

        if (hpW < best_hp){ best_hp=hpW; y_best=R.y_order; }
        // 观测模式：不提前 break，完整记录全轮
    }

    fout << "# SUMMARY best_weighted_HPWL="<<best_hp<<"\n";
    fout.close();
    cout<<"[Observe] wrote metrics to "<<out_path<<endl;
    return 0;
}
