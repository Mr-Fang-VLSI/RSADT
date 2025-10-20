#include "ocCapTDAG.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv){
    if(argc<3){
        std::cout<<"Usage: ./build/oc_captdag m h [dV=1] [T=0->no-limit] [v=2] [K=4]\n";
        return 0;
    }
    int m=std::atoi(argv[1]);
    int h=std::atoi(argv[2]);
    long long dV = (argc>=4)? std::atoll(argv[3]) : 1;
    int T = (argc>=5)? std::atoi(argv[4]) : 0;
    int v = (argc>=6)? std::atoi(argv[5]) : 2;
    int K = (argc>=7)? std::atoi(argv[6]) : 4;

    ocCapTDAG::Config cfg;
    cfg.dV = dV;
    cfg.T = T;
    cfg.vlevel = v;
    cfg.progress = true;
    cfg.logfile = "run.log";
    cfg.log_append = false;
    cfg.max_labels_per_key = K;

    ocCapTDAG solver(cfg);
    std::cout<<"\n=== OC-CapT DAG (APT) m="<<m<<" h="<<h<<" T="<<(T>0?T:m*h)<<" v="<<v<<" K="<<K<<" ===\n";

    try{
        auto R = solver.solve(m,h);
        std::cout<<"Result: HPWL="<<R.hpwl<<", maxΔ="<<R.max_delta<<", OC="<<(R.oc_ok?"OK":"FAIL")<<"\n";
        std::cout<<"(details in run.log)\n";
    }catch(const std::exception& e){
        std::cout<<"ERROR: "<<e.what()<<"\n(details in run.log)\n";
        return 1;
    }
    return 0;
}
