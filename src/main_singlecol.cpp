// src/main_singlecol.cpp
#include "lightOCShortest.h"
#include "momentum_weighter.h"
#include "placement_parser.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <cctype>
#include <limits>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
template<typename T> using Mat = std::vector<std::vector<T>>;

// ---- helpers ----
static bool load_matrix_ll(const string& path, int R, int C, Mat<long long>& M){
    if(path.empty()) return false;
    std::ifstream fin(path);
    if(!fin) return false;
    M.assign(R, std::vector<long long>(C,1));
    long long x; int r=0, c=0;
    while(fin>>x){
        M[r][c]=x; if(++c==C){ c=0; ++r; if(r==R) break; }
    }
    return (r==R && c==0);
}
static bool parse_flag(const string& s, const string& key, string& out){
    auto p=s.find('=');
    if(p==string::npos) return false;
    if(s.substr(0,p)==key){ out=s.substr(p+1); return true; }
    return false;
}
static bool is_integer(const char* s){
    if(!s || !*s) return false;
    // allow optional +/- then digits
    const char* p=s;
    if(*p=='+' || *p=='-') ++p;
    if(!*p) return false;
    while(*p){ if(!std::isdigit((unsigned char)*p)) return false; ++p; }
    return true;
}

int main(int argc,char** argv){
    if(argc<3){
        cout<<"Usage: ./build/oc_singlecol m h [rounds=10] [dV=1]\n"
              "  --T=20               Max-T\n"
              "  --alpha=0.5          momentum alpha\n"
              "  --eta=0.5            momentum step\n"
              "  --top=0.10           top ratio\n"
              "  --cap=32             weight cap\n"
              "  --scale=1000         export scale for integer weights\n"
              "  --wH path            (optional) m x (h-1) horizontal weights\n"
              "  --wV path            (optional) (m-1) x h vertical weights\n"
              "  --pow2=1/0           use 2^k bit-pack (default 1)\n"
              "  --outfile path       output placement file (default auto)\n"
              "  --progress=1/0       verbose\n";
        return 0;
    }

    const int m  = std::atoi(argv[1]);
    const int h  = std::atoi(argv[2]);

    // defaults
    int rounds   = 10;
    long long dV = 1;
    long long T  = std::max(m,h);
    double alpha=0.5, eta=0.5, top=0.10, cap=32.0, scale=1000.0;
    bool progress=true;
    bool use_pow2 = true;
    string wH_path, wV_path, out_path;

    // 1) parse flags first
    for(int i=1;i<argc;++i){
        string s(argv[i]), v;
        if(parse_flag(s,"--T",v))          T = std::llround(std::stod(v));
        else if(parse_flag(s,"--alpha",v)) alpha = std::stod(v);
        else if(parse_flag(s,"--eta",v))   eta   = std::stod(v);
        else if(parse_flag(s,"--top",v))   top   = std::stod(v);
        else if(parse_flag(s,"--cap",v))   cap   = std::stod(v);
        else if(parse_flag(s,"--scale",v)) scale = std::stod(v);
        else if(parse_flag(s,"--rounds",v))rounds= std::stoi(v);
        else if(parse_flag(s,"--dV",v))    dV    = std::llround(std::stod(v));
        else if(parse_flag(s,"--pow2",v))  use_pow2 = (v!="0");
        else if(parse_flag(s,"--progress",v)) progress=(v=="1");
        else if(s=="--wH" && i+1<argc)     wH_path=argv[++i];
        else if(s=="--wV" && i+1<argc)     wV_path=argv[++i];
        else if(s=="--outfile" && i+1<argc) out_path=argv[++i];
    }
    // 2) then apply numeric positional overrides safely
    if(argc>=4 && is_integer(argv[3])) rounds = std::atoi(argv[3]);
    if(argc>=5 && is_integer(argv[4])) dV     = std::atoll(argv[4]);

    if(dV==0){
        cerr<<"[Warn] dV parsed as 0; reset to 1 to avoid degenerate DP.\n";
        dV = 1;
    }

    // weights
    Mat<long long> wH(m, std::vector<long long>(h-1, (long long)std::llround(scale)));
    Mat<long long> wV(m-1, std::vector<long long>(h,   (long long)std::llround(scale)));
    if(!wH_path.empty()){
        if(!load_matrix_ll(wH_path, m, h-1, wH)) cerr<<"[Warn] load wH failed, fallback to 1.\n";
    }
    if(!wV_path.empty()){
        if(!load_matrix_ll(wV_path, m-1, h, wV)) cerr<<"[Warn] load wV failed, fallback to 1.\n";
    }

    MomentumWeighter::Params WP;
    WP.alpha = alpha; WP.eta = eta; WP.top_ratio = top; WP.cap_max = cap; WP.export_scale = scale;
    MomentumWeighter weighter(m,h,WP);
    weighter.set_initial(wH, wV, scale);

    cout<<"[SingleCol-TD] m="<<m<<" h="<<h<<" N="<<(m*h)
        <<"  rounds="<<rounds<<"  dV="<<dV<<"  T="<<T
        <<"  alpha="<<alpha<<" eta="<<eta<<" top="<<top
        <<"  pow2="<<(use_pow2?1:0)<<endl;

    Mat<int> y_best; long long best_hp=std::numeric_limits<long long>::max();

    for(int it=1; it<=rounds; ++it){
        // export integer weights
        Mat<long long> WH_ll, WV_ll;
        weighter.export_ll(WH_ll, WV_ll);

        lightOCShortest::Config cfg;
        cfg.dV = dV;                       // *** critical: avoid 0 ***
        cfg.verbose = progress;
        cfg.progress = progress;
        cfg.max_steps = -1;
        cfg.omp_threads = 0;
        cfg.low_mem = true;
        cfg.use_pow2_pack = use_pow2;
        cfg.W.enabled = true;
        cfg.W.wH = WH_ll;
        cfg.W.wV = WV_ll;

        lightOCShortest solver(cfg);

        auto t0 = std::chrono::high_resolution_clock::now();
        auto R  = solver.solve(m,h);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double,std::milli>(t1-t0).count();

        long long hpW = lightOCShortest::hpwl_sum_weighted(R.y_order, cfg.W, 1);
        long long hpE = lightOCShortest::hpwl_sum_equal   (R.y_order, 1);

        cout<<"  [Round "<<it<<"] solve_ms="<<ms
            <<"  HPWL_w="<<hpW<<"  HPWL_eq="<<hpE
            <<"  report="<<R.total_cost<<endl;

        bool ok = lightOCShortest::check_OC_lin(R.y_order);
        if(!ok) cerr<<"  [OC] FAIL\n";

        // timing-driven update with T
        auto rep = weighter.update_from_layout(R.y_order, T, /*pick_by_ratio=*/true);
        if(progress){
            cout<<"  [Updt] WNS="<<rep.WNS<<"  TNS="<<rep.TNS
                <<"  maxL="<<rep.max_edge_len
                <<"  picked="<<rep.picked_edges<<endl;
        }

        if(hpW < best_hp){ best_hp=hpW; y_best=R.y_order; }

        if(rep.WNS >= 0){
            y_best  = R.y_order;       // <-- 新增
            best_hp = hpW;             //    （可选，不影响正确性）
            std::cout << "  [TD] reach target: all edges <= T." << std::endl;
            break;
        }
    }

    if(y_best.empty()){
        cerr<<"[SingleCol-TD] ERROR: no solution captured.\n";
        return 2;
    }

    // output placement file
    if(out_path.empty()){
        out_path = PlacementParser::default_filename(m,h,T);
    }
    string header = PlacementParser::default_header(m,h,/*cols=*/1);
    bool ok_out = PlacementParser::dump_single_column(
        y_best, m, h, out_path, header, /*x_index=*/0, /*y_from_one=*/true
    );
    if(!ok_out){
        cerr<<"[SingleCol-TD] write output failed.\n";
        return 3;
    }
    cout<<"[SingleCol-TD] wrote placement to: "<< out_path << endl;
    return 0;
}
