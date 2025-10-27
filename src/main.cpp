#include "oc_dp.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <string>
#include <iomanip>

static bool load_matrix(const std::string& path, int R, int C, std::vector<std::vector<long long>>& M){
    if(path.empty()) return false;
    std::ifstream fin(path);
    if(!fin) return false;
    M.assign(R, std::vector<long long>(C, 0));
    long long x; int r=0,c=0;
    while(fin>>x){
        M[r][c]=x;
        if(++c==C){ c=0; ++r; if(r==R) break; }
    }
    return (r==R && c==0);
}

static int getFlagInt(int argc,char**argv,const std::string& key,int defv){
    std::string k="--"+key+"=";
    for(int i=1;i<argc;++i){
        std::string s(argv[i]);
        if(s.rfind(k,0)==0) return std::atoi(s.substr(k.size()).c_str());
    }
    return defv;
}
static bool getFlagBool(int argc,char**argv,const std::string& key,bool defv){
    int v = getFlagInt(argc,argv,key, defv?1:0);
    return v!=0;
}
static std::string getFlagStr(int argc,char**argv,const std::string& key,const std::string& defv){
    std::string k="--"+key;
    for(int i=1;i<argc;++i){
        std::string s(argv[i]);
        if(s==k && i+1<argc) return std::string(argv[i+1]);
        if(s.rfind(k+"=",0)==0) return s.substr(k.size()+1);
    }
    return defv;
}

int main(int argc,char** argv){
    if(argc<3){
        std::cout<<"Usage: "<<argv[0]<<" m h"
                 <<"\n  [--pow2=1] [--progress=1] [--verbose=0] [--small=64]"
                 <<"\n  [--verify=1] [--check-oc=1]"
                 <<"\n  [--wH path] [--wV path] [--print=0]\n";
        return 0;
    }
    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);

    // 配置
    OCSolveConfig cfg;
    cfg.use_pow2_pack = getFlagBool(argc,argv,"pow2", true);
    cfg.progress      = getFlagBool(argc,argv,"progress", true);
    cfg.verbose       = getFlagBool(argc,argv,"verbose", false);
    cfg.small_thresh  = getFlagInt (argc,argv,"small", 64);

    // 权重
    OCWeights W;
    std::string wH_path = getFlagStr(argc,argv,"wH","");
    std::string wV_path = getFlagStr(argc,argv,"wV","");
    if(!wH_path.empty() && !wV_path.empty()){
        if(load_matrix(wH_path, m, h-1, W.wH) && load_matrix(wV_path, m-1, h, W.wV)){
            W.enabled=true;
        }else{
            std::cerr<<"[Warn] failed to load weights; fallback to equal 1\n";
            W.enabled=false;
        }
    }

    // 计时 solve()
    OCDP solver(cfg);
    auto t0 = std::chrono::high_resolution_clock::now();
    OCResult R = solver.solve(m,h,W);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double,std::milli>(t1-t0).count();

    // 输出
    std::cout<<"[OC-DP] m="<<m<<" h="<<h<<" N="<<(m*h)
             <<"  pow2="<<(cfg.use_pow2_pack?1:0)
             <<"  time_ms="<<ms<<"\n";
    // 验证（不计时）
    bool verify = getFlagBool(argc,argv,"verify", true);
    bool ckoc   = getFlagBool(argc,argv,"check-oc", true);
    if(verify){
        long long HP = W.enabled ? OCDP::hpwl_sum_weighted(R.y_order, W)
                                 : OCDP::hpwl_sum_equal(R.y_order);
        std::cout<<"[Verify] HPWL="<<HP<<"  (report="<<R.total_cost<<")  "
                 <<(HP==R.total_cost?"OK":"MISMATCH")<<"\n";
    }
    if(ckoc){
        bool ok = OCDP::check_OC(R.y_order);
        std::cout<<"[OC] "<<(ok?"OK":"FAIL")<<"\n";
    }
    // 可选打印矩阵（小规模）
    if(getFlagBool(argc,argv,"print", false)){
        for(int i=0;i<m;++i){
            for(int j=0;j<h;++j){
                std::cout<<std::setw(4)<<R.y_order[i][j];
            }
            std::cout<<"\n";
        }
    }
    return 0;
}
