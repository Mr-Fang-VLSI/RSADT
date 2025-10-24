#include "oc_closure.h"
#include <iostream>
#include <cstdlib>
#include <string>

int main(int argc,char** argv){
    if(argc<3){
        std::cout << "Usage: ./build/oc_closure m h [dV=1]\n"
                     "  --verify=1/0      compute HPWL check (default=1)\n"
                     "  --check-oc=1/0    check monotone OC (default=1)\n"
                     "  --progress=1/0    print Dinkelbach progress (default=0)\n"
                     "  --dump-y path     write y matrix to text file\n";
        return 0;
    }
    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);
    long long dV = (argc>=4 ? std::atoll(argv[3]) : 1); (void)dV;

    OCClosure::Config cfg;
    cfg.verify = true; cfg.check_oc = true; cfg.progress=false;

    for(int i=1;i<argc;++i){
        std::string s(argv[i]); auto pos=s.find('=');
        if(pos!=std::string::npos){
            auto k=s.substr(0,pos), v=s.substr(pos+1);
            if(k=="--verify")     cfg.verify   = (v=="1");
            else if(k=="--check-oc") cfg.check_oc = (v=="1");
            else if(k=="--progress") cfg.progress = (v=="1");
            else if(k=="--dump-y")   cfg.dump_y_path = v;
        }
    }

    std::cout << "[Closure] m="<<m<<" h="<<h<<" N="<<(m*h)
              << " (closure exact, weighted=0)\n";

    OCClosure solver(cfg);
    auto R = solver.solve(m,h);

    std::cout << "[Verify] HPWL="<<R.hpwl<<"\n";
    std::cout << "[OC] " << (R.oc_ok?"OK":"FAIL") << "\n";
    if(!cfg.dump_y_path.empty()){
        std::cout << "[Dump] y -> " << cfg.dump_y_path << "\n";
    }
    return 0;
}
