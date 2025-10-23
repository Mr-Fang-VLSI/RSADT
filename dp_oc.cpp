#include "dp_oc.h"
#include <unordered_map>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <chrono>
#include <sys/resource.h>   // getrusage

using std::vector;
using std::unordered_map;
using hr_clock = std::chrono::steady_clock;

static inline long long llabsll(long long x){ return x>=0?x:-x; }

// ---------- ?????:wV=wH=1 ----------
void DpOcSolver::build_C_unit(int m,int h, vector<vector<long long>>& C){
    C.assign(m, vector<long long>(h, 0));
    // ???:?? top -1,bottom +1,?? 0
    for(int j=0;j<h;++j){
        C[0][j]      += -1;
        C[m-1][j]    += +1;
    }
    // ???:?? left -1,right +1
    for(int i=0;i<m;++i){
        C[i][0]      += -1;
        C[i][h-1]    += +1;
    }
}

// ---------- Base-(h+1) pow ? ----------
vector<uint64_t> DpOcSolver::powB_list(int m, uint64_t B){
    vector<uint64_t> pw(m,1);
    for(int k=m-2; k>=0; --k){
        __uint128_t t = (__uint128_t)pw[k+1]*(__uint128_t)B;
        if(t > (__uint128_t)std::numeric_limits<uint64_t>::max()){
            throw std::runtime_error("state key overflow: try smaller h or implement big key.");
        }
        pw[k] = (uint64_t)t;
    }
    return pw;
}

// ---------- ?? RSS ----------
double DpOcSolver::get_peak_rss_mb(){
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#if defined(__APPLE__) && defined(__MACH__)
    // mac ? ru_maxrss ???;Linux ?? KB
    return ru.ru_maxrss / (1024.0*1024.0);
#else
    return ru.ru_maxrss / 1024.0; // Linux: KB -> MB
#endif
}

// ---------- ??? ----------
DpResult DpOcSolver::solve(int m, int h){
    const int n = m*h;
    const uint64_t B = (uint64_t)h + 1ull;

    vector<vector<long long>> C;
    build_C_unit(m, h, C);

    // pow ?
    vector<uint64_t> pw = powB_list(m, B);

    // ?? key(?? a[i]=h)
    uint64_t end_key = 0;
    for(int i=0;i<m;++i) end_key += pw[i]*(uint64_t)h;

    auto t0 = hr_clock::now();

    unordered_map<uint64_t, long long> cur;
    cur.reserve(4096);
    cur.emplace(0ull, 0ll);

    // ????
    // parents[level+1][key2] = (prev_key, chosen_row)
    vector< unordered_map<uint64_t, std::pair<uint64_t,int>> > parents;
    if(cfg_.reconstruct){
        parents.resize(n+1);
        for(int L=0; L<=n; ++L) parents[L].reserve(4096);
    }

    // ?? DP
    for(int L=0; L<n; ++L){
        unordered_map<uint64_t, long long> nxt;
        // ??????????(?????? m/2 ?????)
        nxt.reserve(std::max<size_t>(16, cur.size()*(std::min(m, h)>>1)));

        for(const auto& kv : cur){
            uint64_t key = kv.first;
            long long cost = kv.second;

            // ?????? i:a[i]<h ? (i==0 || a[i]+1 <= a[i-1]) ????
            for(int i=0;i<m;++i){
                int ai = digit_at(key, i, pw, B);
                if(ai >= h) continue;
                if(i > 0){
                    int a_prev = digit_at(key, i-1, pw, B);
                    if(ai + 1 > a_prev) continue;
                }
                uint64_t key2 = encode_inc(key, i, pw);
                long long new_cost = cost + C[i][ai]*(long long)(L+1);

                auto it = nxt.find(key2);
                if(it==nxt.end() || new_cost < it->second){
                    nxt[key2] = new_cost;
                    if(cfg_.reconstruct){
                        parents[L+1][key2] = {key, i};
                    }
                }
            }
        }

        cur.swap(nxt);
        if(cfg_.verbose && ((L+1) % std::max(1, n/4) == 0)){
            fprintf(stderr, "[DP] level %d/%d, states=%zu\n", L+1, n, cur.size());
        }
    }

    DpResult R;
    auto it = cur.find(end_key);
    if(it == cur.end()){
        throw std::runtime_error("Terminal state not found (unexpected under OC).");
    }
    R.best_hpwl = it->second;

    // ?? / ????
    auto t1 = hr_clock::now();
    std::chrono::duration<double> dt = t1 - t0;
    R.stats.time_sec = dt.count();
    R.stats.peak_rss_mb = get_peak_rss_mb();
    R.stats.final_states = cur.size();

    // ???? y
    if(cfg_.reconstruct){
        R.y.assign(m, vector<int>(h,0));
        uint64_t k = end_key;
        for(int L=n; L>=1; --L){
            auto pr = parents[L].find(k);
            if(pr == parents[L].end()){
                throw std::runtime_error("Parent link missing. Increase parents[] reserve or bug.");
            }
            uint64_t prev_k = pr->second.first;
            int sel_i = pr->second.second;
            int j0 = digit_at(prev_k, sel_i, pw, B);
            R.y[sel_i][j0] = L;
            k = prev_k;
        }
    }
    return R;
}

// ---------- ?? ----------
long long DpOcSolver::hpwl_neighbors(const vector<vector<int>>& y){
    const int m=(int)y.size(), h=(int)y[0].size();
    long long s=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) s += llabsll((long long)y[i][j+1]-y[i][j]);
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j) s += llabsll((long long)y[i+1][j]-y[i][j]);
    return s;
}
long long DpOcSolver::hpwl_boundary(const vector<vector<int>>& y){
    const int m=(int)y.size(), h=(int)y[0].size();
    long long s=0;
    for(int j=0;j<h;++j) s += llabsll((long long)y[m-1][j]-y[0][j]);
    for(int i=0;i<m;++i) s += llabsll((long long)y[i][h-1]-y[i][0]);
    return s;
}
