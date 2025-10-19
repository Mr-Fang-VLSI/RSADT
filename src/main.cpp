#include "ParamClosureOCPlacer.h"
#include <iostream>
#include <chrono>
#include <cstdlib>

using std::cout;
using std::endl;

static void run_case(int m, int h, long long dV) {
    ParamClosureOCPlacer::Config cfg;
    cfg.dV = dV;
    cfg.verbose = true;

    cout << "\n=== Param-Closure Single-Column Test m=" << m << " h=" << h
         << " (dV=" << dV << ") ===\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    ParamClosureOCPlacer placer(cfg);
    auto R = placer.solve(m, h);
    auto t1 = std::chrono::high_resolution_clock::now();

    cout << "y_order (1.." << m*h << "):\n";
    for (int i = 0; i < m; ++i) {
        cout << "  ";
        for (int j = 0; j < h; ++j)
            cout << R.y_order[i][j] << (j+1==h?'\n':' ');
    }
    cout << "[Check] OC=" << (R.oc_ok? "OK":"FAIL")
         << ", Unique=" << (R.unique_ok? "OK":"FAIL")
         << ", cost=" << R.cost
         << ", HPWL=" << R.actual_hpwl
         << ", Diff=" << (R.actual_hpwl - R.cost) << "\n";

    auto us = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
    cout << "[Time] " << us << " us\n";
}

int main(int argc, char** argv) {
    long long dV = 1;
    run_case(3, 3, dV);
    run_case(4, 4, dV);
    run_case(8, 8, dV);   // 这里应得到 472（一列最优）

    // 自定义： ./build/param_closure m h [dV]
    if (argc == 3 || argc == 4) {
        int m = std::atoi(argv[1]);
        int h = std::atoi(argv[2]);
        long long DV = (argc==4 ? std::atoll(argv[3]) : 1);
        run_case(m, h, DV);
    } else {
        cout << "\n(可选) 自定义: ./build/param_closure m h [dV]\n"
             << "  例如 32×32: ./build/param_closure 32 32 1\n";
    }
    return 0;
}
