#include "oc_closure.h"
#include <algorithm>
#include <queue>
#include <limits>
#include <cassert>
#include <iostream>
#include <unordered_set>

using std::vector; using std::string;

// ---------- 验证/统计 ----------
bool OCClosure::check_OC_lin(const vector<vector<int>>& y){
    int m=(int)y.size(), h=(int)y[0].size();
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) if(y[i][j] > y[i][j+1]) return false;
    for(int j=0;j<h;++j) for(int i=0;i+1<m;++i) if(y[i][j] > y[i+1][j]) return false;
    return true;
}
long long OCClosure::hpwl_sum_equal(const vector<vector<int>>& y, long long dV){
    int m=(int)y.size(), h=(int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) S += llabs((long long)y[i][j+1]-y[i][j]) * dV;
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j) S += llabs((long long)y[i+1][j]-y[i][j]) * dV;
    return S;
}
long long OCClosure::hpwl_sum_weighted(const vector<vector<int>>& y,
                                       const OCWeights& W, long long dV){
    if(!W.enabled) return hpwl_sum_equal(y, dV);
    int m=(int)y.size(), h=(int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j)
        S += llabs((long long)y[i][j+1]-y[i][j]) * W.wH[i][j];
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j)
        S += llabs((long long)y[i+1][j]-y[i][j]) * W.wV[i][j];
    return S * dV;
}

// ---------- φ 构造 ----------
static inline long long H_at(const OCWeights& W,int m,int h,int I,int J){
    if(!W.enabled) return (0<=I && I<m && 0<=J && J<h-1)? 1LL:0LL;
    return (0<=I && I<m && 0<=J && J<h-1)? W.wH[I][J] : 0LL;
}
static inline long long V_at(const OCWeights& W,int m,int h,int I,int J){
    if(!W.enabled) return (0<=I && I<m-1 && 0<=J && J<h)? 1LL:0LL;
    return (0<=I && I<m-1 && 0<=J && J<h)? W.wV[I][J] : 0LL;
}
void OCClosure::build_phi(int m,int h, const OCWeights& W){
    phi_.assign(m, vector<long long>(h, 0));
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        long long up   = V_at(W,m,h,i-1,j);
        long long left = H_at(W,m,h,i, j-1);
        long long down = V_at(W,m,h,i,   j);
        long long right= H_at(W,m,h,i,   j);
        long long v = (up + left) - (down + right);
        phi_[i][j] = v * cfg_.dV;
    }
    // 原始 DAG（向右/向下）—— 用于“向上闭合”约束（选 u 必须选其后继）
    m_=m; h_=h; N_=m*h;
    succ_.assign(N_, {});
    auto inb=[&](int x,int lo,int hi){ return (lo<=x && x<hi); };
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        int u=id(i,j);
        if(inb(i+1,0,m)) succ_[u].push_back(id(i+1,j));
        if(inb(j+1,0,h)) succ_[u].push_back(id(i,j+1));
    }
    if(cfg_.verbose){
        std::cout << "[phi] built, m="<<m_<<" h="<<h_<<" N="<<N_<<"\n";
    }
}

// ---------- 小工具 ----------
void OCClosure::complement_mask(const vector<char>& mask, const vector<char>& S, vector<char>& maskLeft){
    int n=(int)mask.size(); maskLeft.assign(n,0);
    for(int i=0;i<n;++i) maskLeft[i] = (mask[i] && !S[i]) ? 1 : 0;
}
long long OCClosure::sum_w_on_mask(const vector<char>& mask) const {
    long long s=0;
    for(int i=0;i<m_;++i) for(int j=0;j<h_;++j){
        int u=id(i,j); if(mask[u]) s += phi_[i][j];
    }
    return s;
}
int OCClosure::count_on_mask(const vector<char>& mask) const{
    int c=0; for(char b:mask) if(b) ++c; return c;
}
void OCClosure::sinks_upward_closed(const vector<char>& mask, vector<char>& S_sink) const{
    S_sink.assign(N_,0);
    bool any=false;
    for(int u=0; u<N_; ++u){
        if(!mask[u]) continue;
        bool hasSucc=false;
        for(int v: succ_[u]) if(mask[v]){ hasSucc=true; break; }
        if(!hasSucc){ S_sink[u]=1; any=true; }
    }
    if(!any){
        for(int i=m_-1;i>=0;--i) for(int j=h_-1;j>=0;--j){
            int u=id(i,j); if(mask[u]){ S_sink[u]=1; return; }
        }
    }
}
bool OCClosure::is_upward_closed_on_mask(const vector<char>& mask, const vector<char>& S) const{
    for(int u=0; u<N_; ++u){
        if(!mask[u] || !S[u]) continue;
        for(int v: succ_[u]){
            if(mask[v] && !S[v]) return false; // u in S but succ v not in S
        }
    }
    return true;
}

// ---------- Dinic for min-cut ----------
struct Dinic {
    struct E { int v; long long cap; int rev; };
    int n=0, s=0, t=0;
    std::vector<std::vector<E>> g;
    std::vector<int> lvl, it;

    explicit Dinic(int n_): n(n_), g(n_), lvl(n_), it(n_) {}
    void addEdge(int u,int v,long long c){
        if(c<0) c=0;
        E a{v,c,(int)g[v].size()};
        E b{u,0,(int)g[u].size()};
        g[u].push_back(a); g[v].push_back(b);
    }
    bool bfs(){
        std::fill(lvl.begin(), lvl.end(), -1);
        std::queue<int> q; lvl[s]=0; q.push(s);
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(const auto& e:g[u]) if(e.cap>0 && lvl[e.v]<0){
                lvl[e.v]=lvl[u]+1; q.push(e.v);
            }
        }
        return lvl[t]>=0;
    }
    long long dfs(int u,long long f){
        if(u==t || f==0) return f;
        for(int &i=it[u]; i<(int)g[u].size(); ++i){
            auto &e=g[u][i];
            if(e.cap>0 && lvl[e.v]==lvl[u]+1){
                long long got=dfs(e.v, std::min(f, e.cap));
                if(got>0){
                    e.cap -= got;
                    g[e.v][e.rev].cap += got;
                    return got;
                }
            }
        }
        return 0;
    }
    long long maxflow(int S,int T){
        s=S; t=T; long long flow=0;
        const long long INF=(std::numeric_limits<long long>::max)()/4;
        while(bfs()){
            std::fill(it.begin(), it.end(), 0);
            while(true){
                long long got=dfs(s, INF);
                if(got==0) break;
                flow += got;
            }
        }
        return flow;
    }
    void reachable_from_source(std::vector<char>& reach){
        int n=(int)g.size();
        reach.assign(n,0);
        std::queue<int> q; q.push(s); reach[s]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(const auto& e:g[u]) if(e.cap>0 && !reach[e.v]){
                reach[e.v]=1; q.push(e.v);
            }
        }
    }
};

// ---------- 最大权闭包：upward-closed ----------
long long OCClosure::max_weight_upward_closure(const vector<char>& mask,
                                               const vector<long long>& q,
                                               vector<char>& maskS) const{
    const int S = N_, T = N_+1;
    Dinic D(N_+2);

    long long sumPos = 0, sumNeg = 0, sumAbs = 0;
    for(int u=0; u<N_; ++u){
        if(!mask[u]) continue;
        long long qi = q[u];
        if(qi>0){ D.addEdge(S,u,qi); sumPos += qi; sumAbs += qi; }
        else if(qi<0){ D.addEdge(u,T,-qi); sumNeg += -qi; sumAbs += -qi; }
    }
    long long INF = sumAbs + 1; if(INF < 2) INF = 2; // 关键修复：必须大于任意可行割容量

    // upward-closed: 选 u 必须选其后继 v
    for(int u=0; u<N_; ++u){
        if(!mask[u]) continue;
        for(int v: succ_[u]){
            if(mask[v]) D.addEdge(u, v, INF);
        }
    }
    (void)D.maxflow(S,T);

    vector<char> reach;
    D.reachable_from_source(reach);
    maskS.assign(N_,0);
    long long sumQ=0;
    for(int u=0; u<N_; ++u){
        if(mask[u] && reach[u]){ maskS[u]=1; sumQ += q[u]; }
    }

    if(cfg_.verbose){
        int cnt=0; for(int u=0;u<N_;++u) if(maskS[u]) ++cnt;
        std::cout << "[closure] |mask|="<<count_on_mask(mask)
                  << " sumPos="<<sumPos<<" sumNeg="<<sumNeg
                  << " sumAbs="<<sumAbs<<" INF="<<INF
                  << "  => |S|="<<cnt<<" sumQ="<<sumQ
                  << " upClosed="<<(is_upward_closed_on_mask(mask,maskS)?"Y":"N")
                  << "\n";
    }

    return sumQ;
}

// ---------- Dinkelbach：最大密度 upward-closed（带“最小基数”字典序择优） ----------
void OCClosure::find_max_density_upward_closed(const vector<char>& mask,
                                               vector<char>& maskS,
                                               long long& A_sum, int& B_count) const{
    // w_i = phi(i,j)
    long long minw = (std::numeric_limits<long long>::max)();
    for(int i=0;i<m_;++i) for(int j=0;j<h_;++j){
        int u=id(i,j); if(!mask[u]) continue;
        minw = std::min(minw, phi_[i][j]);
    }
    long long A = minw - 1; // λ=A/B with B=1
    long long B = 1;

    vector<long long> q(N_,0);
    vector<char> S(N_,0), S_prev(N_,0);

    int guard=0;
    while(true){
        // q_i = w_i * B - A
        for(int i=0;i<m_;++i) for(int j=0;j<h_;++j){
            int u=id(i,j); if(!mask[u]){ q[u]=0; continue; }
            long long wi = phi_[i][j];
            __int128 tmp = (__int128)wi * (__int128)B - (__int128)A;
            if(tmp > (__int128)std::numeric_limits<long long>::max()) tmp = std::numeric_limits<long long>::max();
            if(tmp < (__int128)std::numeric_limits<long long>::min()) tmp = std::numeric_limits<long long>::min();
            q[u] = (long long)tmp;
        }
        long long sumQ = max_weight_upward_closure(mask, q, S);

        if(cfg_.progress){
            long long As=0; int Bs=0;
            for(int u=0;u<N_;++u) if(mask[u] && S[u]){ ++Bs; int ii=u/h_, jj=u%h_; As+=phi_[ii][jj]; }
            std::cout << "[Dinkelbach] iter="<<guard
                      << " B="<<B<<" A="<<A
                      << " sumQ="<<sumQ
                      << " -> |S|="<<Bs<<" sumW="<<As << std::endl;
        }

        if(sumQ==0){
            break; // 已到 λ*
        }
        // 更新 λ = A(S)/|S|
        long long Anew=0; int Bnew=0;
        for(int u=0;u<N_;++u) if(mask[u] && S[u]){
            ++Bnew; int i=u/h_, j=u%h_; Anew += phi_[i][j];
        }
        if(Bnew==0){
            // 极罕见兜底
            sinks_upward_closed(mask, S);
            Anew=0; Bnew=0;
            for(int u=0;u<N_;++u) if(mask[u] && S[u]){ ++Bnew; int i=u/h_, j=u%h_; Anew+=phi_[i][j]; }
            break;
        }
        if(S == S_prev) break; // 数值反复
        S_prev = S;
        A = Anew; B = Bnew;
        if(++guard > 256) break; // 保险
    }

    // —— 字典序择优：在 λ*=A/B 下，取“最小基数”的最优闭合集 ——
    vector<long long> qstar(N_,0);
    long long U=0;
    for(int i=0;i<m_;++i) for(int j=0;j<h_;++j){
        int u=id(i,j); if(!mask[u]){ qstar[u]=0; continue; }
        long long wi=phi_[i][j];
        __int128 t=(__int128)wi*(__int128)B - (__int128)A;
        if(t > (__int128)std::numeric_limits<long long>::max()) t=std::numeric_limits<long long>::max();
        if(t < (__int128)std::numeric_limits<long long>::min()) t=std::numeric_limits<long long>::min();
        qstar[u]=(long long)t;
        U += (qstar[u]>=0? qstar[u] : -qstar[u]);
    }
    long long K = U + 1; if(K<2) K=2;

    vector<long long> qlex(N_,0);
    for(int u=0;u<N_;++u){
        if(!mask[u]){ qlex[u]=0; continue; }
        __int128 t = (__int128)qstar[u] * (__int128)K - (__int128)1;
        if(t > (__int128)std::numeric_limits<long long>::max()) t=std::numeric_limits<long long>::max();
        if(t < (__int128)std::numeric_limits<long long>::min()) t=std::numeric_limits<long long>::min();
        qlex[u]=(long long)t;
    }
    long long sumQlex = max_weight_upward_closure(mask, qlex, S);

    // —— 安全兜底：若 S==mask 或 S 为空，取 sink 集合 ——
    int Bs=0; for(int u=0;u<N_;++u) if(mask[u] && S[u]) ++Bs;
    bool same=true; for(int u=0;u<N_;++u) if((bool)S[u]!=(bool)mask[u]){ same=false; break; }
    if(same || Bs==0){
        if(cfg_.verbose){
            std::cout << "[lex] same="<<same<<" |S|="<<Bs<<" -> fallback to sinks\n";
        }
        sinks_upward_closed(mask, S);
        Bs=0; for(int u=0;u<N_;++u) if(mask[u] && S[u]) ++Bs;
    }

    // 输出
    maskS = S;
    A_sum = 0; B_count = 0;
    for(int u=0; u<N_; ++u) if(mask[u] && S[u]){
        ++B_count; int i=u/h_, j=u%h_; A_sum += phi_[i][j];
    }
    (void)sumQlex;
}

// ---------- 递归：先排剩余，再排最大密度块（最后段） ----------
std::vector<int> OCClosure::solve_recursive(const vector<char>& mask) const{
    int cnt = count_on_mask(mask);
    if(cnt==0) return {};
    if(cnt==1){
        for(int u=0; u<N_; ++u) if(mask[u]) return {u};
        return {};
    }
    if(cfg_.verbose){
        std::cout << "[rec] enter |mask|="<<cnt<<" sumW="<<sum_w_on_mask(mask) << "\n";
    }

    vector<char> S, Lmask;
    long long As=0; int Bs=0;
    find_max_density_upward_closed(mask, S, As, Bs);

    // 保险：若出现空 S（不应发生），取一个 sink
    if(Bs<=0){
        if(cfg_.verbose) std::cout << "[rec] Bs==0 -> sinks\n";
        sinks_upward_closed(mask, S);
        Bs = count_on_mask(S);
    }
    // 再保：若 S==mask（理论上已避免），仍强拆为 sink 集合
    bool same=true; for(int u=0;u<N_;++u) if((bool)S[u]!=(bool)mask[u]){ same=false; break; }
    if(same){
        if(cfg_.verbose) std::cout << "[rec] S==mask -> sinks fallback\n";
        sinks_upward_closed(mask, S);
        Bs = count_on_mask(S);
    }

    complement_mask(mask, S, Lmask);
    int Lcnt = count_on_mask(Lmask);
    if(cfg_.verbose){
        std::cout << "[rec] choose |S|="<<Bs<<" |Left|="<<Lcnt
                  << " (upClosed="<<(is_upward_closed_on_mask(mask,S)?"Y":"N")<<")\n";
    }

    auto left = solve_recursive(Lmask); // 先排剩余
    auto tail = solve_recursive(S);     // 再排该块（放在最后）

    std::vector<int> out;
    out.reserve(left.size()+tail.size());
    out.insert(out.end(), left.begin(), left.end());
    out.insert(out.end(), tail.begin(), tail.end());

    if((int)out.size()!=cnt){
        std::cerr << "[FATAL] recursion size mismatch: out="<<out.size()<<" vs cnt="<<cnt<<"\n";
        // 打印未覆盖/重复的 id 以便定位
        std::vector<int> mark(N_,0);
        for(int u=0;u<N_;++u) if(mask[u]) mark[u]-=1; // 期望 -1
        for(int x: out) mark[x]+=1; // 实际 +1
        int miss=0, dup=0;
        for(int u=0;u<N_;++u){
            if(!mask[u]) continue;
            if(mark[u]==-1){ ++miss; std::cerr<<"  missing u="<<u<<"\n"; }
            else if(mark[u]>0){ ++dup; std::cerr<<"  duplicated u="<<u<<" count="<<mark[u]<<"\n"; }
        }
        std::cerr << "[FATAL] miss="<<miss<<" dup="<<dup<<"\n";
        std::exit(2);
    }
    return out;
}

// ---------- 顶层求解 ----------
OCClosureResult OCClosure::solve(int m, int h, const OCWeights& W){
    build_phi(m, h, W);
    vector<char> full(N_, 1);

    if(cfg_.verbose){
        std::cout << "[Closure] m="<<m<<" h="<<h<<" N="<<N_
                  << " (closure exact, weighted="<<(W.enabled?1:0)<<")\n";
    }

    auto order = solve_recursive(full); // 从早到晚的 id 序列

    if((int)order.size()!=N_){
        std::cerr << "[FATAL] order.size()!=N_: "<<order.size()<<" vs "<<N_<<"\n";
        std::exit(3);
    }

    // 转 rank 矩阵
    vector<vector<int>> y(m, vector<int>(h, 0));
    for(int t=0; t<N_; ++t){
        int u = order[t];
        if(u<0 || u>=N_){
            std::cerr << "[FATAL] u out of range: u="<<u<<" N="<<N_<<"\n";
            std::exit(4);
        }
        int i=u/h_, j=u%h_;
        y[i][j] = t+1;
    }
    long long HP = W.enabled ? hpwl_sum_weighted(y, W, 1) : hpwl_sum_equal(y, 1);

    return OCClosureResult{m,h,N_, cfg_.dV, std::move(y), HP};
}
