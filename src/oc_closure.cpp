#include "oc_closure.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <limits>

using std::vector;
using std::string;
using ll = long long;

static inline ll llabsll(ll x){ return x>=0?x:-x; }

// ---------- 验证工具 ----------
long long OCClosure::hpwl_sum_equal(const vector<std::vector<int>>& y){
    int m=(int)y.size(), h=(int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i)
        for(int j=0;j+1<h;++j) S += llabsll((long long)y[i][j+1]-y[i][j]);
    for(int i=0;i+1<m;++i)
        for(int j=0;j<h;++j)   S += llabsll((long long)y[i+1][j]-y[i][j]);
    return S;
}
bool OCClosure::check_OC_lin(const vector<std::vector<int>>& y){
    int m=(int)y.size(), h=(int)y[0].size();
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) if(y[i][j]>y[i][j+1]) return false;
    for(int j=0;j<h;++j) for(int i=0;i+1<m;++i) if(y[i][j]>y[i+1][j]) return false;
    return true;
}
void OCClosure::dump_y(const vector<std::vector<int>>& y, const string& path){
    if(path.empty()) return;
    std::ofstream fout(path);
    int m=(int)y.size(), h=(int)y[0].size();
    for(int i=0;i<m;++i){
        for(int j=0;j<h;++j){
            fout<<y[i][j]<<(j+1<h?' ':'\n');
        }
    }
}

// ---------- φ 构造（等权）: φ = (up+left) - (down+right) ----------
void OCClosure::build_phi(int m,int h){
    m_=m; h_=h; N_=m_*h_;
    phi_.assign(m_, vector<ll>(h_, 0));
    auto H_at=[&](int I,int J)->ll{
        return (0<=I && I<m_ && 0<=J && J<h_-1)? 1LL:0LL;
    };
    auto V_at=[&](int I,int J)->ll{
        return (0<=I && I<m_-1 && 0<=J && J<h_)? 1LL:0LL;
    };
    for(int i=0;i<m_;++i) for(int j=0;j<h_;++j){
        ll up   = V_at(i-1,j);
        ll left = H_at(i, j-1);
        ll down = V_at(i, j);
        ll right= H_at(i, j);
        phi_[i][j] = (up + left) - (down + right);
    }
}


// ---------- 最大权 upward-closed 闭包 ----------
long long OCClosure::max_weight_upward_closure(
        const std::vector<long long>& q,
        const std::vector<char>& mask,
        std::vector<char>& S)
{
    // DAG: 每个节点 u 有“右、下”两个后继 v
    const int s = N_, t = N_ + 1;
    Dinic D(N_ + 2);

    // ---- 1. 统计当前掩码节点数 M 与 |φ|max ----
    int M = 0;
    long long phi_max = 0;
    for (int u = 0; u < N_; ++u) {
        if (!mask[u]) continue;
        ++M;
        int i = u / h_, j = u % h_;
        long long v = std::llabs(phi_[i][j]);
        if (v > phi_max) phi_max = v;
    }
    if (phi_max == 0) phi_max = 1;  // 防止全零
    long long INF = 3LL * M * M * phi_max + 1;  // 理论下界推导式

    // ---- 2. 建立 S→u、u→T 边 ----
    for (int u = 0; u < N_; ++u) {
        if (!mask[u]) continue;
        long long w = q[u];
        if (w > 0) D.addEdge(s, u, w);
        else if (w < 0) D.addEdge(u, t, -w);
    }

    // ---- 3. 加闭包约束边：u→succ(u) ----
    auto add_if = [&](int u, int v) {
        if (u < 0 || v < 0 || u >= N_ || v >= N_) return;
        if (mask[u] && mask[v]) D.addEdge(u, v, INF);
    };
    for (int i = 0; i < m_; ++i)
        for (int j = 0; j < h_; ++j) {
            int u = idx(i, j);
            if (!mask[u]) continue;
            if (j + 1 < h_) add_if(u, idx(i, j + 1));  // 右
            if (i + 1 < m_) add_if(u, idx(i + 1, j));  // 下
        }

    // ---- 4. 最小割求最大流 ----
    (void)D.maxflow(s, t);

    // ---- 5. 从源可达集合提取闭合集 ----
    auto reach = D.mincut_src_reachable(s);
    S.assign(N_, 0);
    long long sumQ = 0;
    for (int u = 0; u < N_; ++u)
        if (mask[u] && reach[u]) {
            S[u] = 1;
            sumQ += q[u];
        }

    return sumQ;
}


// ---------- mask 上的 sink 集合（upward-closed 且非空真子集） ----------
void OCClosure::sinks_upward_closed(const vector<char>& mask, vector<char>& S){
    S.assign(N_,0);
    int cnt=0;
    for(int i=0;i<m_;++i){
        for(int j=0;j<h_;++j){
            int u=idx(i,j);
            if(!mask[u]) continue;
            bool hasSucc = false;
            if(j+1<h_ && mask[idx(i,j+1)]) hasSucc = true;
            if(i+1<m_ && mask[idx(i+1,j)]) hasSucc = true;
            if(!hasSucc){ S[u]=1; ++cnt; }
        }
    }
    if(cnt==0){
        // 兜底：任选一个 mask 中的点
        for(int u=0;u<N_;++u) if(mask[u]){ S[u]=1; break; }
    }
}

// ---------- Dinkelbach 最大密度 upward-closed ----------
void OCClosure::find_max_density_upward_closed(const vector<ll>& qA,
                                               const vector<char>& mask,
                                               vector<char>& Sout,
                                               bool progress){
    const int MAX_IT = 256;

    // 用 sinks 作为非空 seed
    vector<char> S; sinks_upward_closed(mask, S);

    ll A=0; int B=0;
    for(int u=0;u<N_;++u) if(mask[u] && S[u]){ A += qA[u]; ++B; }
    double lambda = (B>0? (double)A/B : 0.0);

    vector<ll> qShift(N_,0);

    for(int it=0; it<MAX_IT; ++it){
        // qShift = qA - round(lambda)  （整数容量）
        ll lam_i = (ll)std::llround(lambda);
        for(int u=0;u<N_;++u) qShift[u] = mask[u]? (qA[u] - lam_i) : 0;

        max_weight_upward_closure(qShift, mask, S);

        // 新的 A,B
        ll Anew=0; int Bnew=0;
        for(int u=0;u<N_;++u) if(mask[u] && S[u]){ Anew += qA[u]; ++Bnew; }

        if(progress){
            std::cout<<"[Dinkelbach] it="<<it<<" A="<<Anew<<" B="<<Bnew
                     <<" lambda="<<lambda<<"\n";
        }

        // 空集修复（保证非空）
        if(Bnew==0){
            sinks_upward_closed(mask, S);
            Anew=0; Bnew=0;
            for(int u=0;u<N_;++u) if(mask[u] && S[u]){ Anew+=qA[u]; ++Bnew; }
        }

        double lam_new = (Bnew>0? (double)Anew/Bnew : 0.0);

        // 收敛 / 停滞
        if(Bnew==B && Anew==A) break;
        if(std::fabs(lam_new - lambda) < 1e-9) break;

        A=Anew; B=Bnew; lambda = lam_new;
    }

    // 词典序偏好（大B优先）：qlex = qA*(2B) - A
    if(B>0){
        vector<ll> qlex(N_,0);
        for(int u=0;u<N_;++u){
            if(mask[u]) qlex[u] = qA[u]*(ll)(2*B) - A;
        }
        max_weight_upward_closure(qlex, mask, S);
    }
    Sout = S;
}

// ---------- 递归分解 ----------
void OCClosure::solve_recursive(const vector<char>& mask, vector<int>& order){
    int M = count_on_mask(mask);
    if(M==0) return;
    if(M==1){
        for(int u=0;u<N_;++u) if(mask[u]){ order.push_back(u); return; }
    }

    // qA = φ
    vector<ll> qA(N_,0);
    for(int i=0;i<m_;++i) for(int j=0;j<h_;++j){
        int u=idx(i,j);
        if(mask[u]) qA[u]=phi_[i][j];
    }

    vector<char> S; find_max_density_upward_closed(qA, mask, S, cfg_.progress);

    // 统一兜底：禁止空集/整集
    int cntS=0; for(int u=0;u<N_;++u) if(mask[u] && S[u]) ++cntS;
    if(cntS==0 || cntS==M){
        sinks_upward_closed(mask, S);
        cntS=0; for(int u=0;u<N_;++u) if(mask[u] && S[u]) ++cntS;
    }

    // L = mask - S
    vector<char> Lmask(N_,0);
    for(int u=0;u<N_;++u) if(mask[u] && !S[u]) Lmask[u]=1;

    solve_recursive(Lmask, order);
    solve_recursive(S, order);
}

// ---------- 对外求解 ----------
ClosureResult OCClosure::solve(int m,int h){
    build_phi(m,h);

    vector<char> full(N_,1);
    vector<int> order; order.reserve(N_);
    solve_recursive(full, order);

    if((int)order.size() != N_){
        std::cerr << "[BUG] order.size()="<<order.size()<<" but N_="<<N_<<"\n";
        throw std::runtime_error("closure order incomplete");
    }

    vector<std::vector<int>> y(m_, vector<int>(h_,0));
    for(int t=0; t<N_; ++t){
        int u=order[t];
        int i = u / h_, j = u % h_;
        y[i][j] = t+1;
    }

    long long HP = 0;
    bool ok = true;
    if(cfg_.verify)   HP = hpwl_sum_equal(y);
    if(cfg_.check_oc) ok = check_OC_lin(y);
    if(!cfg_.dump_y_path.empty()) dump_y(y, cfg_.dump_y_path);

    return ClosureResult{m_, h_, N_, 1, std::move(y), HP, ok};
}
