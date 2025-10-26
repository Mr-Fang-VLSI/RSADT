#include "reduced_dp.h"
#include <queue>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <functional>

// ---------- 构造 ----------
ReducedDP::ReducedDP(int m, int h, const ReducedDPConfig& cfg)
: m_(m), h_(h), N_(m*h), cfg_(cfg) {
    build_phi();
}

// ---------- phi（等权） ----------
void ReducedDP::build_phi(){
    phi_.assign(m_, std::vector<int>(h_, 0));
    for(int i=0;i<m_;++i){
        for(int j=0;j<h_;++j){
            int up    = (i>0)?1:0;
            int left  = (j>0)?1:0;
            int down  = (i<m_-1)?1:0;
            int right = (j<h_-1)?1:0;
            // 等权：phi = (上+左) - (下+右)
            phi_[i][j] = (up+left) - (down+right);
        }
    }
}

// ---------- HPWL/OC ----------
long long ReducedDP::hpwl_equal(const std::vector<std::vector<int>>& y){
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long W=0;
    for(int i=0;i<m;++i){
        for(int j=0;j<h;++j){
            if(i+1<m) W += llabs((long long)y[i+1][j]-y[i][j]);
            if(j+1<h) W += llabs((long long)y[i][j+1]-y[i][j]);
        }
    }
    return W;
}

bool ReducedDP::check_oc(const std::vector<std::vector<int>>& y){
    int m = (int)y.size(), h=(int)y[0].size(), N=m*h;
    std::vector<std::pair<int,int>> pos(N+1, {-1,-1});
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        int r = y[i][j];
        if(r<1 || r> N) return false;
        pos[r]={i,j};
    }
    std::vector<std::vector<char>> seen(m, std::vector<char>(h,0));
    for(int r=1;r<=N;++r){
        auto [i,j] = pos[r];
        if(i>0 && !seen[i-1][j]) return false;
        if(j>0 && !seen[i][j-1]) return false;
        seen[i][j]=1;
    }
    return true;
}

// ---------- 初始/终止 ----------
ReducedDP::State ReducedDP::start_state() const{ State s; s.a.assign(m_, 0); return s; }
ReducedDP::State ReducedDP::goal_state()  const{ State s; s.a.assign(m_, (uint16_t)h_); return s; }

// ---------- 规则A：列级串行（单列） ----------
bool ReducedDP::make_macro_col(const State& s, int L0, long long& cost, State& t, EdgeInfo& einfo) const{
    int r = (int)s.a.back(); // 最小列
    if(r>=h_) return false;

    bool exists=false;
    if(s.a[0]==(uint16_t)r) exists=true;
    for(int i=1;i<m_;++i){
        if(s.a[i]==(uint16_t)r && s.a[i-1]>(uint16_t)r){ exists=true; break; }
    }
    if(!exists) return false;

    t = s; cost=0; int L=L0;
    auto push = [&](int i,int j){
        int phi = phi_[i][j];
        cost += (long long)(L+1) * (long long)phi;
        ++L;
        t.a[i] = (uint16_t)(t.a[i]+1);
        if(cfg_.verbose){
            std::cerr<<"  [macro-col] ("<<i<<","<<j<<") phi="<<phi<<"  L->"<<L<<"\n";
        }
    };

    if(t.a[0]==(uint16_t)r) push(0,r); // 顶（多为 -1/0）
    for(int i=1;i<m_;++i){
        if(t.a[i]==(uint16_t)r && t.a[i-1]>(uint16_t)r){
            push(i,r); // 中间链（0 或底 +1）
        }
    }
    einfo.type = EdgeInfo::MACRO_COL; einfo.r = r;
    return true;
}

// ---------- 规则B：左列开行（保守） ----------
bool ReducedDP::make_macro_left(const State& s, int L0, long long& cost, State& t, EdgeInfo& einfo) const{
    if(!cfg_.enable_left_chain) return false;
    // 保守触发：顶已经至少一步（否则 0 有可能插入）
    if(s.a[0]==0) return false;

    // 找最上方连续未开行段（中间行 1..m-2）
    int k=0;
    for(int i=1;i<m_-1;++i){
        if(s.a[i]==0 && s.a[i-1]>0) ++k;
        else break;
    }
    if(k<=0) return false;

    t=s; cost=0; int L=L0;
    auto push = [&](int i,int j){ int phi=phi_[i][j]; cost += (long long)(L+1)*phi; ++L; t.a[i]++; };
    for(int off=0; off<k; ++off){
        int i = 1+off;
        if(!(t.a[i]==0 && t.a[i-1]>0)) return false; // 防御
        push(i,0); // φ=-1
        if(cfg_.verbose){
            std::cerr<<"  [macro-left] ("<<i<<",0) phi=-1  L->"<<L<<"\n";
        }
    }
    einfo.type=EdgeInfo::MACRO_LEFT; einfo.k=k;
    return true;
}

// ---------- 规则C：右列收尾（保守） ----------
bool ReducedDP::make_macro_right(const State& s, int L0, long long& cost, State& t, EdgeInfo& einfo) const{
    if(!cfg_.enable_right_chain) return false;
    if(s.a[0]<(uint16_t)(h_-1)) return false;

    int k=0;
    for(int i=1;i<m_-1;++i){
        if(s.a[i]==(uint16_t)(h_-1) && s.a[i-1]>(uint16_t)(h_-1)) ++k;
        else break;
    }
    if(k<=0) return false;

    t=s; cost=0; int L=L0;
    auto push = [&](int i,int j){ int phi=phi_[i][j]; cost += (long long)(L+1)*phi; ++L; t.a[i]++; };
    for(int off=0; off<k; ++off){
        int i = 1+off;
        if(!(t.a[i]==(uint16_t)(h_-1) && t.a[i-1]>t.a[i])) return false;
        push(i, h_-1); // φ=+1
        if(cfg_.verbose){
            std::cerr<<"  [macro-right] ("<<i<<","<<h_-1<<") phi=+1  L->"<<L<<"\n";
        }
    }
    einfo.type=EdgeInfo::MACRO_RIGHT; einfo.k=k;
    return true;
}

// ---------- 单步兜底 ----------
void ReducedDP::expand_unit_steps(const State& s, int L0,
                                  std::vector<std::pair<State,long long>>& outs,
                                  std::vector<EdgeInfo>& eouts) const{
    outs.clear(); eouts.clear();
    for(int i=0;i<m_;++i){
        if(s.a[i]>= (uint16_t)h_) continue;
        if(!step_feasible(s,i)) continue;
        State t = s; int j = (int)t.a[i];
        int phi = phi_[i][j];
        long long c = (long long)(L0+1) * (long long)phi;
        t.a[i] = (uint16_t)(t.a[i]+1);
        outs.emplace_back(t, c);
        EdgeInfo e; e.type=EdgeInfo::STEP; e.row=i; e.col=j; eouts.push_back(e);
    }
}

// ---------- Canonicalization：最大化“安全宏推进”，但**每次至多吃 1 列** ----------
void ReducedDP::canonicalize(State& s, int L0, long long& delta, std::vector<EdgeInfo>& seq) const{
    delta = 0; seq.clear();
    State cur = s; int L = L0;

    int col_used_total = 0; // 本次 canonicalize 吃了多少列
    bool did_any;

    while(true){
        did_any = false;

        // 1) 左列负相位（保守），可连续
        if(cfg_.enable_left_chain){
            while(true){
                long long c=0; State t; EdgeInfo e;
                if(!make_macro_left(cur, L, c, t, e)) break;
                delta += c; cur = t; L = level_of(cur);
                seq.push_back(e); did_any = true;
            }
        }
        // 2) 右列正相位（保守），可连续
        if(cfg_.enable_right_chain){
            while(true){
                long long c=0; State t; EdgeInfo e;
                if(!make_macro_right(cur, L, c, t, e)) break;
                delta += c; cur = t; L = level_of(cur);
                seq.push_back(e); did_any = true;
            }
        }
        // 3) 列级级联：**每次 canonicalize 最多吃 cfg_.max_col_cascade 列**
        if(cfg_.macro_col_in_canon && col_used_total < cfg_.max_col_cascade){
            long long c=0; State t; EdgeInfo e;
            if(make_macro_col(cur, L, c, t, e)){
                delta += c; cur = t; L = level_of(cur);
                seq.push_back(e); did_any = true;
                ++col_used_total;
            }
        }

        if(!did_any) break;
        // 注意：下一轮循环还会继续尝试左/右保守串，
        // 但“列级级联”总量受 col_used_total 限制（<= max_col_cascade），
        // 因而不会再把多列连续吃空。
    }

    s = cur; // 代表状态
}

// ---------- 回溯展开（递归 lambda 用 std::function） ----------
void ReducedDP::reconstruct_y(const std::vector<EdgeInfo>& path_edges, std::vector<std::vector<int>>& y) const{
    y.assign(m_, std::vector<int>(h_, 0));
    std::vector<uint16_t> a(m_,0);
    int rank=1;

    auto put = [&](int i,int j){ y[i][j]=rank++; a[i]++; };

    std::function<void(const EdgeInfo&)> apply_one = [&](const EdgeInfo& e){
        if(e.type==EdgeInfo::STEP){
            put(e.row, e.col);
        }else if(e.type==EdgeInfo::MACRO_COL){
            int r = e.r;
            if(a[0]==(uint16_t)r) put(0,r);
            for(int i=1;i<m_;++i){
                if(a[i]==(uint16_t)r && a[i-1]>(uint16_t)r) put(i,r);
            }
        }else if(e.type==EdgeInfo::MACRO_LEFT){
            int k = e.k;
            for(int off=0; off<k; ++off){
                int i=1+off;
                if(!(a[0]>0 && a[i]==0 && a[i-1]>0)) continue;
                put(i,0);
            }
        }else if(e.type==EdgeInfo::MACRO_RIGHT){
            int k = e.k;
            for(int off=0; off<k; ++off){
                int i=1+off;
                if(!(a[i]==(uint16_t)(h_-1) && a[i-1]>a[i])) continue;
                put(i, h_-1);
            }
        }else if(e.type==EdgeInfo::COMPOSITE){
            const auto& seq = seq_store_[e.seq_id];
            for(const auto& z: seq) apply_one(z);
        }
    };

    for(const auto& e: path_edges) apply_one(e);
}

void ReducedDP::print_layout(const std::vector<std::vector<int>>& y){
    int m=(int)y.size(), h=(int)y[0].size();
    for(int i=0;i<m;++i){
        for(int j=0;j<h;++j){
            std::cout<<std::setw(4)<<y[i][j];
        }
        std::cout<<"\n";
    }
}

// 组合数近似（用 log Gamma）
static inline double lgamma1p(double x){ return std::lgamma(x+1.0); }
double ReducedDP::binom_approx(int n, int k){
    if(k<0||k>n) return 0.0;
    double lnC = lgamma1p(n) - lgamma1p(k) - lgamma1p(n-k);
    return std::exp(lnC);
}

// ---------- 求解 ----------
ReducedDPResult ReducedDP::solve(){
    using Pair = std::pair<State,long long>;

    State S0 = start_state();
    State ST = goal_state();

    // 分层桶 + best
    std::unordered_map<State,long long,StateHasher> best;
    std::unordered_map<State,EdgeInfo,StateHasher>  last_edge;
    std::unordered_map<State,State,StateHasher>     parent;
    std::vector<std::vector<State>> layers(N_+1);

    size_t state_cnt=0, edge_cnt=0;

    // 起点：谨慎 canonicalization（不会多列连吃）
    State S0c = S0; long long d0=0; std::vector<EdgeInfo> seq0;
    if(cfg_.canonize_start){
        canonicalize(S0c, 0, d0, seq0);
    }
    int L0c = level_of(S0c);
    long long base0 = d0;

    if(L0c>0 && !seq0.empty()){
        EdgeInfo e; e.type=EdgeInfo::COMPOSITE; e.seq_id=(int)seq_store_.size();
        seq_store_.push_back(seq0);
        parent[S0c]=S0; last_edge[S0c]=e; ++edge_cnt;
    }

    best[S0c]=base0;
    layers[L0c].push_back(S0c);
    ++state_cnt;

    for(int L=L0c; L<=N_; ++L){
        auto& cur = layers[L];
        if(cur.empty()) continue;

        // 去重
        std::sort(cur.begin(), cur.end(), [](const State& a, const State& b){ return a.a < b.a; });
        cur.erase(std::unique(cur.begin(), cur.end(), [](const State& a, const State& b){ return a.a==b.a; }), cur.end());

        if(cfg_.verbose)
            std::cerr<<"[Layer "<<L<<"] states="<<cur.size()<<"\n";

        for(const State& s: cur){
            auto it=best.find(s);
            if(it==best.end()) continue;
            long long base = it->second;
            if(s.a==ST.a) continue;

            // 1) 单列宏（产生 1 列的折叠），随后对目标态再做 canonicalize（但受 max_col_cascade 限制）
            {
                long long c=0; State t; EdgeInfo e;
                if(make_macro_col(s, L, c, t, e)){
                    long long d=0; std::vector<EdgeInfo> seq;
                    canonicalize(t, level_of(t), d, seq);
                    long long val = base + c + d;
                    int L2 = level_of(t);

                    std::vector<EdgeInfo> pack; pack.push_back(e);
                    pack.insert(pack.end(), seq.begin(), seq.end());
                    EdgeInfo ec; ec.type=EdgeInfo::COMPOSITE; ec.seq_id=(int)seq_store_.size();
                    seq_store_.push_back(pack);

                    auto jt = best.find(t);
                    if(jt==best.end() || val < jt->second){
                        if(jt==best.end()) { layers[L2].push_back(t); ++state_cnt; }
                        best[t]=val; parent[t]=s; last_edge[t]=ec; ++edge_cnt;
                    }
                }
            }
            // 2) 左/右保守宏
            {
                long long c=0; State t; EdgeInfo e;
                if(make_macro_left(s, L, c, t, e)){
                    long long d=0; std::vector<EdgeInfo> seq;
                    canonicalize(t, level_of(t), d, seq);
                    long long val = base + c + d;
                    int L2 = level_of(t);

                    std::vector<EdgeInfo> pack; pack.push_back(e);
                    pack.insert(pack.end(), seq.begin(), seq.end());
                    EdgeInfo ec; ec.type=EdgeInfo::COMPOSITE; ec.seq_id=(int)seq_store_.size();
                    seq_store_.push_back(pack);

                    auto jt = best.find(t);
                    if(jt==best.end() || val < jt->second){
                        if(jt==best.end()) { layers[L2].push_back(t); ++state_cnt; }
                        best[t]=val; parent[t]=s; last_edge[t]=ec; ++edge_cnt;
                    }
                }
            }
            {
                long long c=0; State t; EdgeInfo e;
                if(make_macro_right(s, L, c, t, e)){
                    long long d=0; std::vector<EdgeInfo> seq;
                    canonicalize(t, level_of(t), d, seq);
                    long long val = base + c + d;
                    int L2 = level_of(t);

                    std::vector<EdgeInfo> pack; pack.push_back(e);
                    pack.insert(pack.end(), seq.begin(), seq.end());
                    EdgeInfo ec; ec.type=EdgeInfo::COMPOSITE; ec.seq_id=(int)seq_store_.size();
                    seq_store_.push_back(pack);

                    auto jt = best.find(t);
                    if(jt==best.end() || val < jt->second){
                        if(jt==best.end()) { layers[L2].push_back(t); ++state_cnt; }
                        best[t]=val; parent[t]=s; last_edge[t]=ec; ++edge_cnt;
                    }
                }
            }
            // 3) 单步兜底 + 立刻 canonicalize（受 max_col_cascade 限制）
            {
                std::vector<Pair> outs; std::vector<EdgeInfo> eouts;
                expand_unit_steps(s, L, outs, eouts);
                for(size_t k=0;k<outs.size();++k){
                    State t = outs[k].first;
                    long long c = outs[k].second;
                    long long d=0; std::vector<EdgeInfo> seq;
                    canonicalize(t, level_of(t), d, seq);
                    long long val = base + c + d;
                    int L2 = level_of(t);

                    std::vector<EdgeInfo> pack; pack.push_back(eouts[k]);
                    pack.insert(pack.end(), seq.begin(), seq.end());
                    EdgeInfo ec; ec.type=EdgeInfo::COMPOSITE; ec.seq_id=(int)seq_store_.size();
                    seq_store_.push_back(pack);

                    auto jt = best.find(t);
                    if(jt==best.end() || val < jt->second){
                        if(jt==best.end()) { layers[L2].push_back(t); ++state_cnt; }
                        best[t]=val; parent[t]=s; last_edge[t]=ec; ++edge_cnt;
                    }
                }
            }
        }
    }

    ReducedDPResult R;
    auto it=best.find(ST);
    if(it==best.end()){
        std::cerr<<"[Error] goal state unreachable.\n";
        return R;
    }
    R.total_cost = it->second;
    R.num_states = state_cnt;
    R.num_edges  = edge_cnt;
    R.states_over_m2h2 = (double)state_cnt / (double)(m_*m_*h_*h_);
    R.states_over_binom = (double)state_cnt / binom_approx(m_+h_, m_);

    // 回溯
    std::vector<EdgeInfo> path;
    State cur = ST;
    while(!(cur.a==S0c.a)){
        EdgeInfo e = last_edge[cur];
        path.push_back(e);
        cur = parent[cur];
    }
    std::reverse(path.begin(), path.end());
    reconstruct_y(path, R.y_order);

    // 打印
    std::cout<<"[ReducedDP+Canon] m="<<m_<<" h="<<h_<<" N="<<N_
             <<"  states="<<R.num_states
             <<"  edges="<<R.num_edges
             <<"  states/(m^2 h^2)="<<std::fixed<<std::setprecision(6)<<R.states_over_m2h2
             <<"  states/C(m+h,m)="<<std::setprecision(6)<<R.states_over_binom
             <<"\n";
    std::cout<<"[ReducedDP+Canon] total_cost="<<R.total_cost<<"\n";

    long long HP = hpwl_equal(R.y_order);
    std::cout<<"[Verify] HPWL="<<HP<<"  (report="<<R.total_cost<<")  "<<(HP==R.total_cost?"OK":"MISMATCH")<<"\n";

    bool ocok = check_oc(R.y_order);
    std::cout<<"[OC] "<<(ocok?"OK":"FAIL")<<"\n";

    if(cfg_.print_layout){
        std::cout<<"[Layout]\n";
        print_layout(R.y_order);
    }
    return R;
}
