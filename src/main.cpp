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
