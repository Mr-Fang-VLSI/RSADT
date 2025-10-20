#include "greedyNoT.h"
#include <fstream>
#include <limits>
#include <algorithm>

namespace {
std::ofstream g_log;
int g_vl = 2;
inline void LOGOPEN2(const std::string& path, bool append, int v){
    if(!g_log.is_open()){
        g_log.open(path, append?(std::ios::out|std::ios::app):(std::ios::out|std::ios::trunc));
        g_vl=v;
    }
}
template<typename F> inline void LOG2(int lvl, F&& f){
    if(g_log.is_open() && g_vl>=lvl){ f(g_log); g_log.flush(); }
}
}

GreedyNoTResult GreedyNoTSolver::solve(int m, int h){
    LOGOPEN2(cfg_.logfile, cfg_.log_append, cfg_.vlevel);
    const int n = m*h;
    LOG2(1, [&](std::ofstream& o){
        o<<"[Greedy-NoT] m="<<m<<" h="<<h<<" n="<<n<<" v="<<cfg_.vlevel<<"\n";
    });

    std::vector<std::vector<int>> y(m, std::vector<int>(h,0));
    std::vector<int> a(m, 0); // skyline：每行已放置的列数（候选为 (i, a[i])）
    long long lin_cost = 0;

    for(int L=1; L<=n; ++L){
        int best_i = -1;
        long long best_w = std::numeric_limits<long long>::max();

        // 候选：i==0 或 a[i] < a[i-1]，且 a[i] < h
        for(int i=0;i<m;++i){
            int j = a[i];
            if(j>=h) continue;
            if(i>0 && !(a[i] < a[i-1])) continue;
            long long w = weight_ij(i, j, m, h);
            // 取 w 最小；打平按 i 从小到大稳定
            if(w < best_w || (w==best_w && i < best_i)){
                best_w = w;
                best_i = i;
            }
        }
        if(best_i<0){
            LOG2(1, [&](std::ofstream& o){ o<<"[Greedy-NoT] no feasible candidate at L="<<L<<"\n"; });
            throw std::runtime_error("Greedy-NoT: no feasible candidate (should not happen for OC)");
        }

        int i = best_i;
        int j = a[i];
        y[i][j] = L;
        lin_cost += (long long)L * weight_ij(i,j,m,h) * cfg_.dV;
        a[i]++;

        if(cfg_.vlevel>=3 && (L<=16 || L==n)){
            LOG2(3, [&](std::ofstream& o){
                o<<"[Greedy-NoT] L="<<L<<" pick ("<<i<<","<<j<<") w="<<weight_ij(i,j,m,h)<<"\n";
            });
        }
    }

    long long hpwl = hpwl_neighbors(y, cfg_.dV);
    bool ocok = check_OC(y);

    LOG2(1, [&](std::ofstream& o){
        o<<"[Greedy-NoT] HPWL="<<hpwl<<" linear="<<lin_cost
         <<" eq="<<(hpwl==lin_cost?"YES":"NO")
         <<" OC="<<(ocok?"OK":"FAIL")<<"\n";
    });
    if(cfg_.vlevel>=3){
        LOG2(3, [&](std::ofstream& o){
            o<<"[y_order]\n";
            for(int i=0;i<m;++i){
                o<<"  ";
                for(int j=0;j<h;++j){
                    o<<(j?" ":"")<<(y[i][j]<10?"  ":" ")<<y[i][j];
                }
                o<<"\n";
            }
        });
    }

    GreedyNoTResult R;
    R.m=m; R.h=h; R.n=n; R.dV=cfg_.dV;
    R.y=std::move(y);
    R.hpwl=hpwl;
    R.linear_cost=lin_cost;
    R.oc_ok=ocok;
    return R;
}
