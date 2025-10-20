#include "ocCapTDAG.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv){
    if(argc < 5){
        std::cerr << "Usage: " << argv[0] << " m h dV T [v=2] [K=0]\n";
        return 1;
    }
    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);
    long long dV = std::atoll(argv[3]);
    int T = std::atoi(argv[4]);
    int v = (argc>=6? std::atoi(argv[5]) : 2);
    int K = (argc>=7? std::atoi(argv[6]) : 0);

    ocCapTDAG::Config cfg;
    cfg.dV = dV;
    cfg.T = T;
    cfg.vlevel = v;
    cfg.logfile = "run.log";
    cfg.log_append = false;
    cfg.max_labels_per_key = K;
    cfg.progress = true;

    try{
        ocCapTDAG solver(cfg);
        auto R = solver.solve(m,h);
        std::cout << "\n=== OC-CapT DAG (APT) m="<<m<<" h="<<h<<" T="<<T
                  <<" v="<<v<<" K="<<K<<" ===\n";
        std::cout << "Result: HPWL="<<R.hpwl
                  <<", maxΔ="<<R.max_delta
                  <<", OC="<<(R.oc_ok?"OK":"FAIL")<<"\n"
                  <<"(details in run.log)\n";
        return 0;
    }catch(const std::exception& e){
        std::cout << "\n=== OC-CapT DAG (APT) m="<<m<<" h="<<h<<" T="<<T
                  <<" v="<<v<<" K="<<K<<" ===\n";
        std::cout << "ERROR: " << e.what() << "\n(details in run.log)\n";
        return 2;
    }
}
