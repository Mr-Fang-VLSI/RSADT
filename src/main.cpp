#include "ParamThresholdRefineOCPlacer.h"
#include <iostream>
#include <chrono>
#include <cstdlib>

using std::cout;
using std::endl;

static void run_case(int m, int h, long long dV, long long SCALE=1000) {
    ParamThresholdRefineOCPlacer::Config cfg;
    cfg.dV = dV;
    cfg.verbose = true;
    cfg.LAM_SCALE = SCALE;

    cout << "\n=== λ-Refine (Max-Closure Chain) Single-Column Test m=" << m
         << " h=" << h << " (dV=" << dV << ", SCALE=" << SCALE << ") ===\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    ParamThresholdRefineOCPlacer placer(cfg);
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
    run_case(8, 8, dV);   // 目标 472

    // 自定义： ./build/lam_refine m h [dV] [SCALE]
    if (argc >= 3) {
        int m = std::atoi(argv[1]);
        int h = std::atoi(argv[2]);
        long long DV = (argc>=4 ? std::atoll(argv[3]) : 1);
        long long SCALE = (argc>=5 ? std::atoll(argv[4]) : 1000);
        run_case(m, h, DV, SCALE);
    } else {
        cout << "\n(可选) 自定义: ./build/lam_refine m h [dV] [SCALE]\n"
             << "  例如 32×20: ./build/lam_refine 32 20 1 1000\n";
    }
    return 0;
}
