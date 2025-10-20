#include "lightOCMaxT.h"
#include "networkOCMaxT.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char** argv){
    if(argc < 3){
        std::cout << "Usage: ./build/oc_maxt m h [dV=1] [T=m*h] [threads=0] [progress=0/1] [mode=2/3] [vlevel=1..3]\n";
        return 0;
    }
    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);
    long long dV = (argc>=4)? std::atoll(argv[3]) : 1;
    int T = (argc>=5)? std::atoi(argv[4]) : m*h;
    int threads = (argc>=6)? std::atoi(argv[5]) : 0;
    bool progress = (argc>=7)? (std::atoi(argv[6])!=0) : true;
    int mode = (argc>=8)? std::atoi(argv[7]) : 2;
    int vlevel = (argc>=9)? std::atoi(argv[8]) : 1;

    std::cout << "\n=== OC-MaxT Test m="<<m<<" h="<<h<<" T="<<T
              << " thr="<<threads<<" mode="<<mode<<" v="<<vlevel<<" ===\n";

    try{
        if(mode==3){
            networkOCMaxT::Config c2;
            c2.dV = dV;
            c2.progress = progress;
            c2.verbose_level = vlevel;
            c2.logfile = "run.log";
            c2.log_append = false;
            // 建议：给每个 akey 一个温和的前沿上限以稳住内存（8x8 可设 2048/4096）
            c2.cap_per_bucket   = 0;   // 如需更稳设 2048
            c2.enable_dominance = true;
            c2.log_every_levels = 1;

            networkOCMaxT solver2(c2);
            auto R = solver2.solve(m,h,T);
            std::cout << "[mode=3] Result: HPWL="<<R.total_cost<<", maxΔ≤T, OC=OK\n";
        }else{
            lightOCMaxT::Config cfg;
            cfg.dV = dV;
            cfg.progress = progress;
            cfg.prefix_strategy = mode; // 兼容你的旧实现
            cfg.omp_threads = threads;
            cfg.verbose_level = vlevel;
            cfg.logfile = "run.log";
            cfg.log_append = false;
            lightOCMaxT solver(cfg);
            auto R = solver.solve(m,h,T);
            std::cout << "[mode!=3] Result: HPWL="<<R.total_cost<<", maxΔ≤T, OC=OK\n";
        }
    }catch(const std::exception& e){
        std::cout << "ERROR: " << e.what() << "\n";
        std::cout << "(See run.log for details)\n";
        return 1;
    }

    std::cout << "(Details logged to run.log)\n";
    return 0;
}
