#include "lightOCShortest.h"
#include <iostream>
#include <cstdlib>
#include <set>
using std::cout; using std::endl;
// ======= 调试核验辅助 =======

static long long hpwl_neighbors(const std::vector<std::vector<int>>& y, long long dV=1) {
    int m=y.size(), h=y[0].size();
    long long S=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) S += std::llabs((long long)y[i][j+1]-y[i][j])*dV;
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j)   S += std::llabs((long long)y[i+1][j]-y[i][j])*dV;
    return S;
}
// RSAD Eq.(2) 边界口径：只看上下行与左右列；四角抵消
static long long hpwl_boundary(const std::vector<std::vector<int>>& y, long long dV=1){
    int m=y.size(), h=y[0].size();
    long long S=0;
    // 顶/底两行的横向差分
    for(int j=0;j+1<h;++j){
        S += std::llabs((long long)y[m-1][j+1]-y[m-1][j])*dV; // top row: i=m-1
        S += std::llabs((long long)y[0][j+1]  -y[0][j])  *dV; // bottom row: i=0
    }
    // 左/右两列的纵向差分
    for(int i=0;i+1<m;++i){
        S += std::llabs((long long)y[i+1][0]  -y[i][0])  *dV; // left col:  j=0
        S += std::llabs((long long)y[i+1][h-1]-y[i][h-1])*dV; // right col: j=h-1
    }
    return S;
}
// 正确的边界化简（端点差），无绝对值：OC保证非负
static long long hpwl_boundary_endpoints(const std::vector<std::vector<int>>& y, long long dV=1){
    int m=y.size(), h=y[0].size();
    long long S=0;
    for(int j=0;j<h;++j)   S += (long long)y[m-1][j] - (long long)y[0][j];   // 顶行 - 底行
    for(int i=0;i<m;++i)   S += (long long)y[i][h-1] - (long long)y[i][0];   // 右列 - 左列
    return S * dV;
}

static void dump_y(const std::vector<std::vector<int>>& y){
    int m=y.size(), h=y[0].size();
    for(int i=m-1;i>=0;--i){ // 从顶行往下打，便于和论文图一致
        for(int j=0;j<h;++j){
            printf("%4d", y[i][j]);
        }
        printf("\n");
    }
}
static void sanity_check(const std::vector<std::vector<int>>& y){
    int m=y.size(), h=y[0].size(), n=m*h;
    std::set<int> s; s.clear();
    int mn=1e9, mx=-1e9, zeros=0;
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        int v=y[i][j];
        if(v==0) ++zeros;
        else { s.insert(v); mn=std::min(mn,v); mx=std::max(mx,v); }
    }
    printf("[Sanity] unique_nonzero=%zu, min=%d, max=%d, zeros=%d\n",
           s.size(), (s.empty()?0:mn), (s.empty()?0:mx), zeros);
}

static void run_case(int m,int h,long long dV,int max_steps){
    lightOCShortest::Config cfg;
    cfg.dV = dV;
    cfg.verbose = true;
    cfg.progress = true;
    cfg.max_steps = max_steps;

    cout << "\n=== OC-Shortest Test m="<<m<<" h="<<h<<" T="<<(max_steps>0?max_steps:m*h)<<" ===\n";
    lightOCShortest solver(cfg);
    auto R = solver.solve(m,h);
    cout << "Result: cost="<<R.total_cost
         << ", HPWL="<<R.hpwl
         << ", OC="<<(R.oc_ok?"OK":"FAIL") << endl;
    // ======= 打印与三口径验算 =======
printf("\n[y_order top->bottom]\n"); 
dump_y(R.y_order);
sanity_check(R.y_order);

long long Hn = hpwl_neighbors(R.y_order, dV);
long long Hb = hpwl_boundary_endpoints(R.y_order, dV);
printf("[HPWL check] neighbors=%lld, boundary=%lld, solver_cost=%lld\n", 
       Hn, Hb, R.total_cost /* 或 R.mcmf_cost，看你的结构 */);
if(Hn!=Hb){
    printf("!! Mismatch between neighbor-sum and boundary-form. Need investigate (indexing/OC/metric).\n");
}

}

int main(int argc,char** argv){
    long long dV = 1;
    run_case(4,4,dV,-1);

    if(argc==3 || argc==4 || argc==5){
        int m = std::atoi(argv[1]);
        int h = std::atoi(argv[2]);
        long long DV = (argc>=4 ? std::atoll(argv[3]) : 1);
        int T = (argc==5 ? std::atoi(argv[4]) : -1); // Max-T 前缀
        run_case(m,h,DV,T);
    }else{
        cout << "\nUsage: ./build/oc_shortest m h [dV] [T]\n";
    }
    return 0;
}
