#include "lightOCShortest.h"
#include <iostream>
#include <cstdlib>

static void run_case(int m, int h, long long dV){
    lightOCShortest::Config cfg;
    cfg.dV = dV;
    cfg.verbose = true;

    std::cout << "\n=== OC-Shortest Test m="<<m<<" h="<<h<<" ===\n";
    lightOCShortest solver(cfg);
    auto R = solver.solve(m,h);

    std::cout << "Result: cost=" << R.total_cost
              << ", HPWL=" << R.hpwl
              << ", OC=" << (R.oc_ok?"OK":"FAIL") << "\n";
}

int main(int argc, char** argv){
    long long dV = 1;
    // sanity
    run_case(4,4,dV);

    if(argc==3 || argc==4){
        int m = std::atoi(argv[1]);
        int h = std::atoi(argv[2]);
        long long DV = (argc==4? std::atoll(argv[3]) : 1);
        run_case(m,h,DV);
    } else {
        std::cout << "\nUsage: ./build/oc_shortest m h [dV]\n";
    }
    return 0;
}
