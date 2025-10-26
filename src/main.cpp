#include "reduced_dp.h"
#include <iostream>
#include <cstdlib>
#include <string>

static int getFlagInt(int argc,char**argv,const std::string& key,int defv){
    std::string k="--"+key+"=";
    for(int i=1;i<argc;++i){
        std::string s(argv[i]);
        if(s.rfind(k,0)==0){
            return std::atoi(s.substr(k.size()).c_str());
        }
    }
    return defv;
}
static bool getFlagBool(int argc,char**argv,const std::string& key,bool defv){
    int v=getFlagInt(argc,argv,key, defv?1:0);
    return v!=0;
}

int main(int argc,char** argv){
    if(argc<3){
        std::cout<<"Usage: "<<argv[0]<<" m h"
                 <<" [--verbose=0/1] [--print_layout=1/0]"
                 <<" [--left_chain=1/0] [--right_chain=1/0]"
                 <<" [--canon_start=1/0]"
                 <<" [--macro_col_canon=1/0] [--max_col_cascade=1]\n";
        return 0;
    }
    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);

    ReducedDPConfig cfg;
    cfg.verbose            = getFlagBool(argc,argv,"verbose", false);
    cfg.print_layout       = getFlagBool(argc,argv,"print_layout", true);
    cfg.enable_left_chain  = getFlagBool(argc,argv,"left_chain", true);
    cfg.enable_right_chain = getFlagBool(argc,argv,"right_chain", true);
    cfg.canonize_start     = getFlagBool(argc,argv,"canon_start", true);
    cfg.macro_col_in_canon = getFlagBool(argc,argv,"macro_col_canon", true);
    cfg.max_col_cascade    = getFlagInt(argc,argv,"max_col_cascade", 1);

    ReducedDP solver(m,h,cfg);
    auto R = solver.solve();
    return 0;
}
