#include "dp_oc.h"
#include <cstdio>
#include <cstring>
#include <iostream>

static void print_y(const std::vector<std::vector<int>>& y){
    for(const auto& row: y){
        for(size_t j=0;j<row.size();++j){
            std::printf("%3d%s", row[j], (j+1==row.size() ? "\n" : " "));
        }
    }
}

int main(int argc, char** argv){
    if(argc < 3){
        std::fprintf(stderr, "Usage: %s <m> <h> [--reconstruct] [--verbose]\n", argv[0]);
        return 1;
    }
    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);
    DpOcSolver::Config cfg;
    for(int i=3;i<argc;++i){
        if(std::strcmp(argv[i],"--reconstruct")==0) cfg.reconstruct=true;
        else if(std::strcmp(argv[i],"--verbose")==0) cfg.verbose=true;
    }

    DpOcSolver solver(cfg);
    auto R = solver.solve(m,h);

    std::cout << "=== DP OC (no Max-T, unweighted) on " << m << "x" << h << " ===\n";
    std::cout << "Best HPWL = " << R.best_hpwl << "\n";
    std::cout << "Time = " << R.stats.time_sec << " s"
              << ", PeakRSS = " << R.stats.peak_rss_mb << " MB"
              << ", FinalStates = " << R.stats.final_states << "\n";

    if(cfg.reconstruct){
        auto hp_nei = DpOcSolver::hpwl_neighbors(R.y);
        auto hp_bnd = DpOcSolver::hpwl_boundary(R.y);
        std::cout << "[Check] HPWL(nei)=" << hp_nei
                  << ", HPWL(boundary)=" << hp_bnd << "\n";
        if(m<=10 && h<=10){
            std::cout << "[y]\n";
            print_y(R.y);
        }
    }
    return 0;
}
