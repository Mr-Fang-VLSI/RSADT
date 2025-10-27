// src/main_td.cpp
#include "lightOCShortest.h"
#include "momentum_weighter.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <cctype>
#include <string>

using std::cout;
using std::cerr;
using std::endl;
using std::string;
template<typename T> using Mat = std::vector<std::vector<T>>;

// ========== ???? ==========
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
static void save_matrix_ll(const string& path, const Mat<long long>& M){
    std::ofstream fout(path);
    for(size_t i=0;i<M.size();++i){
        for(size_t j=0;j<M[i].size();++j){
            fout<<M[i][j]<<(j+1==M[i].size()?'\n':' ');
        }
    }
}

// ????/?? HPWL(?? lightOCShortest ???????;??????????)
static long long hpwl_equal(const Mat<int>& y){
    const int m=(int)y.size(), h=(int)y[0].size();
    long long sum=0;
    for(int i=0;i<m;++i)
        for(int j=0;j<h-1;++j)
            sum += std::llabs((long long)y[i][j+1]- (long long)y[i][j]);
    for(int i=0;i<m-1;++i)
        for(int j=0;j<h;++j)
            sum += std::llabs((long long)y[i+1][j]- (long long)y[i][j]);
    return sum;
}
static bool parse_flag(const string& s, const string& key, string& out){
    auto p=s.find('=');
    if(p==string::npos) return false;
    if(s.substr(0,p)==key){ out=s.substr(p+1); return true; }
    return false;
}
static int argi_of(int argc,char** argv,const string& name){
    for(int i=1;i<argc;++i) if(string(argv[i])==name) return i;
    return -1;
}

// ========== ??? ==========
int main(int argc,char** argv){
    if(argc<3){
        cout<<"Usage: ./build/oc_td m h [rounds=10] [dV=1]\n"
              "  --T=15               Max T ??(????)\n"
              "  --alpha=0.5          ???? a ? [0,1]\n"
              "  --eta=0.5            ?? ? > 0(?? ?=1-a)\n"
              "  --top=0.10           ???? Top ??(? slack ??)\n"
              "  --cap=32             ??????(<=0 ?????)\n"
              "  --scale=1000         ????????\n"
              "  --rounds=10          ????\n"
              "  --wH path            ??????(m x (h-1)),??=1\n"
              "  --wV path            ??????((m-1) x h),??=1\n"
              "  --verify=1/0         ?????? HPWL ??(??=1)\n"
              "  --check-oc=1/0       ?? OC(??=1)\n"
              "  --progress=1/0       ????(??=1)\n";
        return 0;
    }

    const int m  = std::atoi(argv[1]);
    const int h  = std::atoi(argv[2]);
    int rounds   = (argc>=4 ? std::atoi(argv[3]) : 10);
    long long dV = (argc>=5 ? std::atoll(argv[4]) : 1);

    // ????
    long long T = std::max(m,h);       // ?? T
    double alpha=0.5, eta=0.5, top=0.10, cap=32.0, scale=1000.0;
    bool verify=true, check_oc=true, progress=true;
    string wH_path, wV_path;

    for(int i=1;i<argc;++i){
        string s(argv[i]), v;
        if(parse_flag(s,"--T",v)) T = std::llround(std::stod(v));
        else if(parse_flag(s,"--alpha",v)) alpha = std::stod(v);
        else if(parse_flag(s,"--eta",v))   eta   = std::stod(v);
        else if(parse_flag(s,"--top",v))   top   = std::stod(v);
        else if(parse_flag(s,"--cap",v))   cap   = std::stod(v);
        else if(parse_flag(s,"--scale",v)) scale = std::stod(v);
        else if(parse_flag(s,"--rounds",v)) rounds = std::stoi(v);
        else if(parse_flag(s,"--verify",v)) verify=(v=="1");
        else if(parse_flag(s,"--check-oc",v)) check_oc=(v=="1");
        else if(parse_flag(s,"--progress",v)) progress=(v=="1");
        else if(s=="--wH" && i+1<argc) wH_path=argv[++i];
        else if(s=="--wV" && i+1<argc) wV_path=argv[++i];
    }

    // —— ????? —— //
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

    // —— ?? —— //
    cout<<"[TD] m="<<m<<" h="<<h<<" N="<<(m*h)<<"  rounds="<<rounds
        <<"  T="<<T<<"  alpha="<<alpha<<" eta="<<eta<<" top="<<top<<endl;

    Mat<int> y_best; long long best_hp=std::numeric_limits<long long>::max();

    for(int it=1; it<=rounds; ++it){
        // ????????
        Mat<long long> WH_ll, WV_ll;
        weighter.export_ll(WH_ll, WV_ll);

        lightOCShortest::Config cfg;
        cfg.dV = dV;
        cfg.verbose = progress;
        cfg.progress = progress;
        cfg.max_steps = -1;
        cfg.omp_threads = 0;
        cfg.low_mem = true;
        cfg.use_pow2_pack = true;
        cfg.W.enabled = true;
        cfg.W.wH = WH_ll;
        cfg.W.wV = WV_ll;

        lightOCShortest solver(cfg);

        // ?? solve ??
        auto t0 = std::chrono::high_resolution_clock::now();
        auto R  = solver.solve(m,h);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double,std::milli>(t1-t0).count();

        // ?? & ??
        long long hpW = lightOCShortest::hpwl_sum_weighted(R.y_order, cfg.W, 1);
        long long hpE = lightOCShortest::hpwl_sum_equal   (R.y_order, 1);

        if(verify){
            cout<<"[Round "<<it<<"] solve_ms="<<ms
                <<"  HPWL_w="<<hpW<<"  HPWL_eq="<<hpE
                <<"  report="<<R.total_cost<<endl;
        }

        if(check_oc){
            bool ok = lightOCShortest::check_OC_lin(R.y_order);
            if(!ok) cerr<<"  [OC] FAIL\n";
        }

        // ??? layout ?????
        auto rep = weighter.update_from_layout(R.y_order, T, /*pick_by_ratio=*/true);

        if(progress){
            cout<<"  [Updt] WNS="<<rep.WNS<<"  TNS="<<rep.TNS
                <<"  maxL="<<rep.max_edge_len
                <<"  picked="<<rep.picked_edges<<endl;
        }

        // ????
        if(hpW < best_hp){ best_hp=hpW; y_best=R.y_order; }

        // ?? T ????????
        if(rep.WNS >= 0){
            cout<<"  [TD] reach target: all edges <= T."<<endl;
            break;
        }
    }

    // —— ???? —— //
    if(!y_best.empty()){
        cout<<"[TD] best_weighted_HPWL="<<best_hp<<endl;
        if(check_oc){
            bool ok = lightOCShortest::check_OC_lin(y_best);
            cout<<"[OC] "<<(ok?"OK":"FAIL")<<endl;
        }
    }
    return 0;
}
