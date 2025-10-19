#include "KClosureExactPlacer.h"
#include <iostream>
#include <chrono>
#include <cstdlib>
using std::cout; using std::endl;

static void run_case(int m,int h,long long dV){
    KClosureExactPlacer::Config cfg; cfg.dV=dV; cfg.verbose=true;
    cout<<"\n=== K-Closure-Exact Test m="<<m<<" h="<<h<<" (dV="<<dV<<") ===\n";
    auto t0=std::chrono::high_resolution_clock::now();
    KClosureExactPlacer solver(cfg);
    auto R = solver.solve(m,h);
    auto t1=std::chrono::high_resolution_clock::now();

    cout<<"y_order (1.."<<m*h<<"):\n";
    for(int i=0;i<m;++i){ cout<<"  ";
        for(int j=0;j<h;++j) cout<<R.y_order[i][j]<<(j+1==h?'\n':' ');
    }
    cout<<"[Check] OC="<<(R.oc_ok?"OK":"FAIL")
        <<", Unique="<<(R.unique_ok?"OK":"FAIL")
        <<", cost="<<R.cost
        <<", HPWL="<<R.actual_hpwl
        <<", Diff="<<(R.actual_hpwl-R.cost)<<"\n";
    auto ms=std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
    cout<<"[Time] "<<ms<<" ms\n";
}

int main(int argc,char** argv){
    run_case(3,3,1);   // 预期 24
    run_case(4,4,1);   // 预期 60
    run_case(8,8,1);   // 预期 472

    if(argc==3 || argc==4){
        int m=std::atoi(argv[1]); int h=std::atoi(argv[2]);
        long long dV=(argc==4? std::atoll(argv[3]):1);
        run_case(m,h,dV);
    }else{
        cout<<"\n(可选) 自定义: ./build/kclosure_exact m h [dV]\n";
    }
    return 0;
}
