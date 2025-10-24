#pragma once
#include <vector>
#include <string>
#include <cstdint>

struct ClosureResult {
    int m=0, h=0, n=0;
    long long dV=1;
    std::vector<std::vector<int>> y_order; // rank in [1..n]
    long long hpwl = 0;
    bool oc_ok = false;
};

class OCClosure {
public:
    struct Config {
        long long dV = 1;         // 等权版本未使用，但保留接口
        bool verify = true;
        bool check_oc = true;
        bool progress = false;
        std::string dump_y_path;
    };
    explicit OCClosure(const Config& cfg): cfg_(cfg) {}

    ClosureResult solve(int m, int h);

private:
    Config cfg_;

    // φ (m x h), 等权HPWL差分: (up+left) - (down+right)
    std::vector<std::vector<long long>> phi_;
    void build_phi(int m,int h);

    // 验证/工具
    static long long hpwl_sum_equal(const std::vector<std::vector<int>>& y);
    static bool check_OC_lin(const std::vector<std::vector<int>>& y);
    static void dump_y(const std::vector<std::vector<int>>& y, const std::string& path);

    // ------- 闭包求解（最大权闭包 & 最大密度 upward-closed） -------
    struct Dinic {
        struct Edge { int v; long long cap; int rev; };
        int N;
        std::vector<std::vector<Edge>> G;
        std::vector<int> level, it;
        Dinic(int n=0):N(n),G(n),level(n),it(n){}
        void reset(int n){ N=n; G.assign(n,{}); level.assign(n,0); it.assign(n,0); }
        void addEdge(int u,int v,long long c){
            Edge a{v,c,(int)G[v].size()};
            Edge b{u,0,(int)G[u].size()};
            G[u].push_back(a); G[v].push_back(b);
        }
        bool bfs(int s,int t){
            std::fill(level.begin(),level.end(),-1);
            std::vector<int> q; q.reserve(N);
            level[s]=0; q.push_back(s);
            for(size_t i=0;i<q.size();++i){
                int u=q[i];
                for(const auto& e:G[u]) if(e.cap>0 && level[e.v]<0){
                    level[e.v]=level[u]+1; q.push_back(e.v);
                }
            }
            return level[t]>=0;
        }
        long long dfs(int u,int t,long long f){
            if(u==t) return f;
            for(int &i=it[u]; i<(int)G[u].size(); ++i){
                auto &e=G[u][i];
                if(e.cap>0 && level[e.v]==level[u]+1){
                    long long ret=dfs(e.v,t,std::min(f,e.cap));
                    if(ret>0){
                        e.cap -= ret;
                        G[e.v][e.rev].cap += ret;
                        return ret;
                    }
                }
            }
            return 0;
        }
        long long maxflow(int s,int t){
            long long flow=0, aug;
            while(bfs(s,t)){
                std::fill(it.begin(),it.end(),0);
                while((aug=dfs(s,t,(std::numeric_limits<long long>::max)()))>0) flow+=aug;
            }
            return flow;
        }
        std::vector<char> mincut_src_reachable(int s){
            std::vector<char> vis(N,0);
            std::vector<int> st; st.push_back(s); vis[s]=1;
            while(!st.empty()){
                int u=st.back(); st.pop_back();
                for(const auto& e:G[u]) if(e.cap>0 && !vis[e.v]){
                    vis[e.v]=1; st.push_back(e.v);
                }
            }
            return vis;
        }
    };

    int m_=0, h_=0, N_=0;
    inline int idx(int i,int j) const { return i*h_ + j; }

    // 最大权 upward-closed 闭包
    long long max_weight_upward_closure(const std::vector<long long>& q,
                                    const std::vector<char>& mask,
                                    std::vector<char>& S);


    // mask 上的 sink 集合（upward-closed 且非空真子集）
    void sinks_upward_closed(const std::vector<char>& mask, std::vector<char>& S);

    int count_on_mask(const std::vector<char>& mask) const {
        int c=0; for(int u=0;u<N_;++u) if(mask[u]) ++c; return c;
    }

    // Dinkelbach：在 upward-closed 集类上最大化 sum(qA)/|S|
    void find_max_density_upward_closed(const std::vector<long long>& qA,
                                        const std::vector<char>& mask,
                                        std::vector<char>& Sout,
                                        bool progress);

    // 递归分解：mask -> left + S
    void solve_recursive(const std::vector<char>& mask, std::vector<int>& order);
};
