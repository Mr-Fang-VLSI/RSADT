#include "lightOCMaxT.h"
#include <iostream>
#include <cstdlib>
using std::cout; using std::endl;

static void run_case(int m,int h,long long dV,int T,int threads,bool progress){
    lightOCMaxT::Config cfg;
    cfg.dV = dV;
    cfg.verbose = true;
    cfg.progress = progress;
    cfg.omp_threads = threads;
    cfg.low_mem = true;

    cout << "\n=== OC-MaxT Test m="<<m<<" h="<<h
         << " T="<<(T>0?T:m*h)
         << " thr="<<(threads>0?threads:0) <<" ===\n";
    lightOCMaxT solver(cfg);
    auto R = solver.solve(m,h,T);

    if (R.T == R.n_full) {
        cout << "Result: cost="<<R.total_cost
             << ", HPWL="<<R.hpwl_full
             << ", OC="<<(R.oc_ok?"OK":"FAIL") << endl;
    } else {
        cout << "Result: prefix_cost(T="<<R.T<<")="<<R.total_cost
             << ", HPWL_prefix="<<R.hpwl_prefix
             << ", OC="<<(R.oc_ok?"OK":"FAIL") << endl;
    }
}

int main(int argc,char** argv){
    // Usage: ./build/oc_maxt m h [dV] [T] [threads] [progress:0/1]
    if(argc>=3){
        int m = std::atoi(argv[1]);
        int h = std::atoi(argv[2]);
        long long DV = (argc>=4 ? std::atoll(argv[3]) : 1);
        int T  = (argc>=5 ? std::atoi(argv[4]) : -1);  // -1 -> full
        int thr= (argc>=6 ? std::atoi(argv[5]) : 0);
        bool prog = (argc>=7 ? std::atoi(argv[6])!=0 : true);
        run_case(m,h,DV,T,thr,prog);
    } else {
        cout << "Usage: ./build/oc_maxt m h [dV] [T] [threads] [progress:0/1]\n";
        // 简单 sanity
        run_case(4,4,1,-1,0,true);
        run_case(8,20,1,80,0,true);
    }
    return 0;
}
