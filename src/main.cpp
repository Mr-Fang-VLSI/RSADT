#include "oc_closure.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdlib>
using namespace std;

static bool load_matrix(const string& path, int R, int C, vector<vector<long long>>& M){
    if(path.empty()) return false;
    ifstream fin(path);
    if(!fin) return false;
    M.assign(R, vector<long long>(C,0));
    long long x; int r=0, c=0;
    while(fin>>x){
        M[r][c]=x; if(++c==C){ c=0; ++r; if(r==R) break; }
    }
    return (r==R && c==0);
}

int main(int argc,char** argv){
    if(argc<3){
        cout << "Usage: ./build/oc_closure m h [dV=1]\n"
                "  --verify=1/0      compute HPWL check (default=1)\n"
                "  --check-oc=1/0    check monotone OC (default=1)\n"
                "  --progress=1/0    print Dinkelbach progress (default=0)\n"
                "  --verbose=1/0     print detailed recursion/flow logs (default=1)\n"
                "  --wH path         horizontal weights file: m x (h-1)\n"
                "  --wV path         vertical   weights file: (m-1) x h\n"
                "  --dump-y path     write y matrix to text file\n";
        return 0;
    }
    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);
    long long dV = (argc>=4 ? std::atoll(argv[3]) : 1);

    bool verify=true, check_oc=true, progress=false, verbose=true;
    string wH_path, wV_path, dump_path;

    for(int i=1;i<argc;++i){
        string s(argv[i]); auto pos=s.find('=');
        auto eat=[&](const string& k,const string& v){
            if(k=="--verify") verify=(v=="1");
            else if(k=="--check-oc") check_oc=(v=="1");
            else if(k=="--progress") progress=(v=="1");
            else if(k=="--verbose") verbose=(v=="1");
        };
        if(pos!=string::npos) eat(s.substr(0,pos), s.substr(pos+1));
        else if(s=="--wH" && i+1<argc) wH_path=argv[++i];
        else if(s=="--wV" && i+1<argc) wV_path=argv[++i];
        else if(s=="--dump-y" && i+1<argc) dump_path=argv[++i];
    }

    OCWeights W;
    if(!wH_path.empty() && !wV_path.empty()){
        if(load_matrix(wH_path, m, h-1, W.wH) && load_matrix(wV_path, m-1, h, W.wV)){
            W.enabled = true;
        }else{
            cerr<<"[Warn] failed to load weights; fallback to equal 1\n";
            W.enabled = false;
        }
    }else{
        W.enabled = false;
    }

    OCClosure::Config cfg; cfg.dV=dV; cfg.verbose=verbose; cfg.progress=progress;
    OCClosure solver(cfg);

    auto t0 = std::chrono::high_resolution_clock::now();
    auto R  = solver.solve(m,h,W);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double,std::milli>(t1-t0).count();

    cout << "Solve-only time: " << ms << " ms\n";

    if(verify){
        long long HP = W.enabled ? OCClosure::hpwl_sum_weighted(R.y_order, W, 1)
                                 : OCClosure::hpwl_sum_equal(R.y_order, 1);
        cout << "[Verify] HPWL="<<HP<<" (solver_report="<<R.total_cost<<")\n";
    }
    if(check_oc){
        bool ok = OCClosure::check_OC_lin(R.y_order);
        cout << "[OC] " << (ok?"OK":"FAIL") << "\n";
    }
    if(!dump_path.empty()){
        ofstream fo(dump_path);
        for(int i=0;i<m;++i){
            for(int j=0;j<h;++j){
                fo << R.y_order[i][j] << (j+1==h?'\n':' ');
            }
        }
        fo.close();
        cout << "[Dump] y -> " << dump_path << "\n";
    }
    return 0;
}
