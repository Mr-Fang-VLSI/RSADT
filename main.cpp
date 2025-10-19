#include "IdaStarColumnPlacer.h"
#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;

static void run_case(int m, int h, long long dV) {
    IdaStarColumnPlacer::Config cfg;
    cfg.dV = dV;
    cfg.verbose = true;

    cout << "\n=== IDA* Single-Column Test m=" << m << " h=" << h
         << " (dV=" << dV << ") ===\n";
    IdaStarColumnPlacer placer(cfg);
    auto R = placer.solve(m, h);

    // pretty print y_order
    cout << "y_order (ranks):\n";
    for (int i = m-1; i >= 0; --i) {
        for (int j = 0; j < h; ++j) {
            cout << R.y_order[i][j] << (j+1==h?'\n':'\t');
        }
    }
}

int main(int argc, char** argv) {
    long long dV = 1;
    // demo
    run_case(4, 4, dV);
    // run_case(8, 8, dV);

    // optional CLI: ./build/ida_star m h [dV]
    if (argc == 3 || argc == 4) {
        int m = std::atoi(argv[1]);
        int h = std::atoi(argv[2]);
        long long DV = (argc==4 ? std::atoll(argv[3]) : 1);
        run_case(m, h, DV);
    } else {
        cout << "\n(Optional) custom: ./build/ida_star m h [dV]\n";
    }
    return 0;
}
