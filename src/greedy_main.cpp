#include "greedyNoT.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv){
    if(argc<3){
        std::cout<<"Usage: ./build/oc_greedy m h [dV=1] [v=2]\n";
        return 0;
    }
    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);
    long long dV = (argc>=4)? std::atoll(argv[3]) : 1;
    int v = (argc>=5)? std::atoi(argv[4]) : 2;

    std::cout<<"\n=== Greedy-NoT test m="<<m<<" h="<<h<<" dV="<<dV<<" v="<<v<<" ===\n";
    GreedyNoTSolver::Config cfg;
    cfg.dV = dV;
    cfg.vlevel = v;
    cfg.logfile = "run.log";
    cfg.log_append = false;

    try{
        GreedyNoTSolver solver(cfg);
        auto R = solver.solve(m,h);
        std::cout<<"Result: HPWL="<<R.hpwl<<" linear="<<R.linear_cost
                 <<" eq="<<(R.hpwl==R.linear_cost?"YES":"NO")
                 <<" OC="<<(R.oc_ok?"OK":"FAIL")<<"\n";
        std::cout<<"(Details in run.log)\n";
    }catch(const std::exception& e){
        std::cout<<"ERROR: "<<e.what()<<"\n(See run.log)\n";
        return 1;
    }
    return 0;
}
