#include "multicol_splitter.h"
#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>

using std::string;

static bool parse_flag(const string& s, const string& key, string& out){
    auto p = s.find('=');
    if (p == string::npos) return false;
    if (s.substr(0,p) == key){
        out = s.substr(p+1);
        return true;
    }
    return false;
}

int main(int argc, char** argv){
    if (argc < 4){
        std::cout
            << "Usage: ./build/oc_multicol M H rounds "
            << "[--K=2|4] [--T=12] [--policy=soft|lagrange|irl1] "
            << "[--alpha=0.9 --eta=1.0 --gamma=0.3 --beta=3.0] "
            << "[--rho=0.3] [--lam=2.0 --delta=2.0] "
            << "[--freeze=1] [--cap=64 --scale=1000] [--dV=1] [--pow2=1] [--progress=1] "
            << "[--homotopy=1 --sigma=0.8 --shrink=0.9] [--snake=0|1] "
            << "[--single_bin=build/oc_policy_singlecol] [--work_dir=.]"
            << std::endl;
        return 0;
    }

    MultiColOptions opt;
    opt.m       = std::atoi(argv[1]);
    opt.h       = std::atoi(argv[2]);
    opt.rounds  = std::atoi(argv[3]);

    string v;
    for (int i=4;i<argc;++i){
        string s(argv[i]);
        if      (parse_flag(s,"--K",v))          opt.K          = std::stoi(v);
        else if (parse_flag(s,"--T",v))          opt.T          = std::stoi(v);
        else if (parse_flag(s,"--policy",v))     opt.policy     = v;
        else if (parse_flag(s,"--alpha",v))      opt.alpha      = std::stod(v);
        else if (parse_flag(s,"--eta",v))        opt.eta        = std::stod(v);
        else if (parse_flag(s,"--gamma",v))      opt.gamma      = std::stod(v);
        else if (parse_flag(s,"--beta",v))       opt.beta       = std::stod(v);
        else if (parse_flag(s,"--rho",v))        opt.rho        = std::stod(v);
        else if (parse_flag(s,"--lam",v))        opt.lam        = std::stod(v);
        else if (parse_flag(s,"--delta",v))      opt.delta      = std::stod(v);
        else if (parse_flag(s,"--freeze",v))     opt.freeze     = std::stoi(v);
        else if (parse_flag(s,"--cap",v))        opt.cap        = std::stoll(v);
        else if (parse_flag(s,"--scale",v))      opt.scale      = std::stoll(v);
        else if (parse_flag(s,"--dV",v))         opt.dV         = std::stoi(v);
        else if (parse_flag(s,"--pow2",v))       opt.pow2       = std::stoi(v);
        else if (parse_flag(s,"--progress",v))   opt.progress   = std::stoi(v);
        else if (parse_flag(s,"--homotopy",v))   opt.homotopy   = std::stoi(v);
        else if (parse_flag(s,"--sigma",v))      opt.sigma      = std::stod(v);
        else if (parse_flag(s,"--shrink",v))     opt.shrink     = std::stod(v);
        else if (parse_flag(s,"--snake",v))      opt.snake      = std::stoi(v);
        else if (parse_flag(s,"--single_bin",v)) opt.single_bin = v;
        else if (parse_flag(s,"--work_dir",v))   opt.work_dir   = v;
    }

    std::cout << "[MultiCol] m="<<opt.m<<" h="<<opt.h
              <<" K="<<opt.K<<" T="<<opt.T
              <<" rounds="<<opt.rounds
              <<" policy="<<opt.policy
              <<" snake="<<opt.snake
              <<" work_dir="<<opt.work_dir << std::endl;

    MultiColSplitter splitter;
    auto t0 = std::chrono::high_resolution_clock::now();
    std::string out_path = splitter.run(opt);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double,std::milli>(t1-t0).count();

    if (out_path.empty()){
        std::cerr << "[MultiCol] FAILED to generate placement." << std::endl;
        return 1;
    }
    std::cout << "[MultiCol] done. placement="<<out_path
              <<"  total_ms="<<ms << " ms" << std::endl;
    return 0;
}
