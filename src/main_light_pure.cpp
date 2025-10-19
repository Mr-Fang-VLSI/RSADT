// #include "PureMcmfColumnPlacer.h"
// #include <iostream>
// #include <cstdlib>

// using std::cout;
// using std::endl;

// static void run_case(int m, int h, long long dV) {
//     PureMcmfColumnPlacer::Config cfg;
//     cfg.dV = dV;
//     cfg.verbose = true;

//     cout << "\n=== Pure-MCMF Single-Column Test m=" << m << " h=" << h
//          << " (dV=" << dV << ") ===\n";
//     PureMcmfColumnPlacer placer(cfg);
//     auto R = placer.solve(m, h);

//     cout << "x_of_col: ";
//     for (int j = 0; j < h; ++j) cout << R.x_of_col[j] << (j+1==h?'\n':' ');
//     cout << "y_order (global labels 1.." << (m*h) << "):\n";
//     for (int i = 0; i < m; ++i) {
//         cout << "  ";
//         for (int j = 0; j < h; ++j) cout << R.y_order[i][j] << (j+1==h?'\n':' ');
//     }
//     cout << "[Check] OC=" << (R.oc_ok ? "OK" : "FAIL")
//          << ", Unique=" << (R.unique_ok ? "OK" : "FAIL")
//          << ", MCMF cost=" << R.mcmf_cost
//          << ", Actual HPWL=" << R.actual_hpwl
//          << ", Diff=" << (R.actual_hpwl - R.mcmf_cost) << "\n";
// }

// int main(int argc, char** argv) {
//     // 你的验证：dH 巨大 ⇒ 全部进一根 site 列；这里直接不使用 dH，x 全 0。
//     long long dV = 1;
//     // run_case(3, 3, dV);
//     run_case(4, 4, dV);
//     // run_case(8, 8, dV);

//     // 可选自定义： ./build/pure_mcmf m h [dV]
//     if (argc == 3 || argc == 4) {
//         int m = std::atoi(argv[1]);
//         int h = std::atoi(argv[2]);
//         long long DV = (argc==4 ? std::atoll(argv[3]) : 1);
//         run_case(m, h, DV);
//     } else {
//         cout << "\n(可选) 自定义: ./build/pure_mcmf m h [dV]\n";
//     }
//     return 0;
// }

#include "lightPureMcmf.h"
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;

static void run_case(int m, int h, long long dV) {
    lightPureMcmf::Config cfg;
    cfg.dV = dV;
    cfg.verbose = true;

    cout << "\n=== LightMCMF Test m=" << m << " h=" << h << " ===\n";
    lightPureMcmf placer(cfg);
    auto R = placer.solve(m, h);
    cout << "MCMF cost=" << R.mcmf_cost
         << ", HPWL=" << R.actual_hpwl
         << ", OC=" << (R.oc_ok?"OK":"FAIL")
         << ", Unique=" << (R.unique_ok?"OK":"FAIL") << endl;
}

int main(int argc, char** argv) {
    long long dV = 1;

    // 默认跑一个小例子，确认可执行
    run_case(4, 4, dV);

    // 自定义： ./build/light_mcmf m h [dV]
    if (argc == 3 || argc == 4) {
        int m = std::atoi(argv[1]);
        int h = std::atoi(argv[2]);
        long long DV = (argc==4 ? std::atoll(argv[3]) : 1);
        run_case(m, h, DV);
    } else {
        cout << "\nUsage: ./build/light_mcmf m h [dV]\n";
    }
    return 0;
}

