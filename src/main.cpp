#include "lightOCShortest.h"
#include <iostream>
#include <cstdlib>

using std::cout; using std::endl;

static void run_case(int m,int h,long long dV,int max_steps){
    lightOCShortest::Config cfg;
    cfg.dV = dV;
    cfg.verbose = true;
    cfg.progress = true;
    cfg.max_steps = max_steps;

    cout << "\n=== OC-Shortest Test m="<<m<<" h="<<h<<" T="<<(max_steps>0?max_steps:m*h)<<" ===\n";
    lightOCShortest solver(cfg);
    auto R = solver.solve(m,h);
    cout << "Result: cost="<<R.total_cost
         << ", HPWL="<<R.hpwl
         << ", OC="<<(R.oc_ok?"OK":"FAIL") << endl;
}

int main(int argc,char** argv){
    long long dV = 1;
    run_case(4,4,dV,-1);

    if(argc==3 || argc==4 || argc==5){
        int m = std::atoi(argv[1]);
        int h = std::atoi(argv[2]);
        long long DV = (argc>=4 ? std::atoll(argv[3]) : 1);
        int T = (argc==5 ? std::atoi(argv[4]) : -1); // Max-T 前缀
        run_case(m,h,DV,T);
    }else{
        cout << "\nUsage: ./build/oc_shortest m h [dV] [T]\n";
    }
    return 0;
}
