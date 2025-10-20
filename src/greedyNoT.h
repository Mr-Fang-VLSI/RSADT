#pragma once
#include <vector>
#include <string>

struct GreedyNoTResult {
    int m=0, h=0, n=0;
    long long dV=1;
    std::vector<std::vector<int>> y;   // 1..n 的次序矩阵
    long long hpwl=0;                  // 邻接 HPWL（按 dV 缩放）
    long long linear_cost=0;           // Σ L·w 的线性代价（应与 hpwl 相等）
    bool oc_ok=true;                   // 偏序检查
};

class GreedyNoTSolver {
public:
    struct Config {
        long long dV = 1;
        int  vlevel = 2;               // 0~3
        std::string logfile = "run.log";
        bool log_append = false;
    };
    explicit GreedyNoTSolver(const Config& c): cfg_(c) {}
    GreedyNoTResult solve(int m, int h);

private:
    static long long weight_ij(int i,int j,int m,int h){
        long long w=0;
        if(i==0)     w -= 1;
        if(i==m-1)   w += 1;
        if(j==0)     w -= 1;
        if(j==h-1)   w += 1;
        return w;
    }
    static long long hpwl_neighbors(const std::vector<std::vector<int>>& y, long long dV){
        const int m=(int)y.size(), h=(int)y[0].size();
        long long S=0;
        auto absl = [](long long x){ return x>=0?x:-x; };
        for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) S += absl((long long)y[i][j+1]-y[i][j])*dV;
        for(int i=0;i+1<m;++i) for(int j=0;j<h;++j)   S += absl((long long)y[i+1][j]-y[i][j])*dV;
        return S;
    }
    static bool check_OC(const std::vector<std::vector<int>>& y){
        const int m=(int)y.size(), h=(int)y[0].size();
        for(int i1=0;i1<m;++i1) for(int j1=0;j1<h;++j1){
            for(int i2=i1;i2<m;++i2) for(int j2=j1;j2<h;++j2)
                if(y[i1][j1] > y[i2][j2]) return false;
        }
        return true;
    }

    Config cfg_;
};
