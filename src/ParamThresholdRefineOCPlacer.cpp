#include "ParamThresholdRefineOCPlacer.h"
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <unordered_set>
#include <queue>
#include <iostream>
#include <algorithm>

using std::vector;
using std::cout;
using std::endl;

static inline int vid(int i, int j, int h) { return i*h + j; }

long long ParamThresholdRefineOCPlacer::hpwl_edges_sum(
    const vector<vector<int>>& y,
    const vector<long long>& x,
    long long dV) {
    (void)x; // 单列 x=0
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long S = 0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j + 1 < h; ++j)
            S += (long long)std::llabs(y[i][j+1] - y[i][j]) * dV; // 水平邻边（dx=0）
    for (int i = 0; i + 1 < m; ++i)
        for (int j = 0; j < h; ++j)
            S += (long long)std::llabs(y[i+1][j] - y[i][j]) * dV; // 垂直邻边
    return S;
}

bool ParamThresholdRefineOCPlacer::check_global_OC(const vector<vector<int>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
    for (int j1 = 0; j1 < h; ++j1)
    for (int i2 = i1; i2 < m; ++i2)
    for (int j2 = j1; j2 < h; ++j2)
        if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

// 单次闭包（min-cut）：容量均为整数，强制 keep 在源侧，block 在汇侧
static int max_closure_with_keep_block(
    int m, int h,
    const vector<long long>& w_scaled, // = base_weight * SCALE
    long long lam,                     // 整数 λ（同尺度）
    const vector<char>& keep,          // 1: 必须在 S(λ)
    const vector<char>& block,         // 1: 必须不在 S(λ)
    vector<char>& inS,                 // out: 源侧闭包
    long long& cuts_run,
    bool verbose)
{
    const int n = m*h;
    const long long INF = 1e12;

    lemon::ListDigraph g;
    auto s = g.addNode();
    auto t = g.addNode();

    vector<lemon::ListDigraph::Node> V(n);
    for (int v = 0; v < n; ++v) V[v] = g.addNode();

    lemon::ListDigraph::ArcMap<long long> cap(g);
    auto add_arc = [&](lemon::ListDigraph::Node u, lemon::ListDigraph::Node v, long long c){
        auto a = g.addArc(u, v); cap[a] = c; return a;
    };

    long long edge_cnt = 0, sumPos = 0;
    for (int v = 0; v < n; ++v) {
        long long p = w_scaled[v] - lam; // 可能为负
        if (p >= 0) { add_arc(s, V[v], p); sumPos += p; ++edge_cnt; }
        else        { add_arc(V[v], t, -p); ++edge_cnt; }
        if (keep[v])  { add_arc(s, V[v], INF); ++edge_cnt; }
        if (block[v]) { add_arc(V[v], t, INF); ++edge_cnt; }
    }
    // 覆盖弧（OC）
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < h; ++j) {
            int v = vid(i,j,h);
            if (i > 0) add_arc(V[v], V[vid(i-1,j,h)], INF), ++edge_cnt;
            if (j > 0) add_arc(V[v], V[vid(i,j-1,h)], INF), ++edge_cnt;
        }

    // TODO(maxT 单列)：若 cfg.maxT_steps>0，可在此加额外 INF 弧约束间距

    if (verbose) {
        cout << "[closure] lam=" << lam
             << " nodes=" << (n+2)
             << " arcs~=" << edge_cnt
             << " keep=" << std::count(keep.begin(), keep.end(), (char)1)
             << " block="<< std::count(block.begin(), block.end(), (char)1) << endl;
    }

    lemon::Preflow<lemon::ListDigraph, lemon::ListDigraph::ArcMap<long long>> pf(g, cap, s, t);
    pf.runMinCut(); ++cuts_run;

    inS.assign(n, 0);
    int sz = 0;
    for (int v = 0; v < n; ++v) {
        bool onS = pf.minCut(V[v]);
        inS[v] = onS ? 1 : 0;
        if (onS) ++sz;
    }
    if (verbose) cout << "  [cut] |S|=" << sz << " sumPos=" << sumPos << endl;
    return sz;
}

PTRResult ParamThresholdRefineOCPlacer::solve(int m, int h) {
    const int n = m*h;
    const long long SCALE = cfg_.LAM_SCALE; // λ 精度
    // w_scaled = base_weight * SCALE
    vector<long long> w_scaled(n, 0);
    for (int i = 0, v=0; i < m; ++i)
        for (int j = 0; j < h; ++j, ++v)
            w_scaled[v] = (long long)base_weight(i,j,m,h) * SCALE;

    // λ 搜索边界（足够宽）：使 S(λ_hi)=∅, S(λ_lo)=V
    const long long WMAX = 2 * SCALE; // base_weight ∈ [-2,2]
    long long lam_hi =  + (WMAX + 10*SCALE);
    long long lam_lo =  - (WMAX + 10*SCALE);

    // 嵌套链端点
    vector<char> keep0(n,0), block0(n,0), Shi, Slo;
    long long cuts = 0;
    int s_hi = max_closure_with_keep_block(m,h,w_scaled,lam_hi,keep0,block0,Shi,cuts,cfg_.verbose); // 预期 0
    int s_lo = max_closure_with_keep_block(m,h,w_scaled,lam_lo,keep0,block0,Slo,cuts,cfg_.verbose); // 预期 n
    (void)s_hi; (void)s_lo;

    // 结果 rank
    vector<int> rank(n, 0); int next_rank = 1;

    // 已确定（赋过 rank）的点，后续强制 keep（安全起见）
    vector<char> fixed(n,0);

    // 递归细化：把 (lam_hi,Shi) 与 (lam_lo,Slo) 区间分到“单点差集”
    std::function<void(long long,const vector<char>&, long long,const vector<char>&)> refine =
    [&](long long lhi, const vector<char>& Shi_, long long llo, const vector<char>& Slo_) {
        // 差集 D = Slo \ Shi
        vector<int> D; D.reserve(n);
        for (int v=0; v<n; ++v) if (!Shi_[v] && Slo_[v] && !fixed[v]) D.push_back(v);
        if (D.empty()) return;
        if ((int)D.size() == 1) {
            int v = D[0];
            rank[v] = next_rank++; fixed[v]=1;
            return;
        }

        // 中点 λ
        long long lmid = lhi + (llo - lhi)/2;
        // 约束：keep = Shi_ U fixed； block = complement(Slo_)
        vector<char> keep(n,0), block(n,0), Smid;
        for (int v=0; v<n; ++v) keep[v] = (Shi_[v] || fixed[v]) ? 1 : 0;
        for (int v=0; v<n; ++v) block[v] = (Slo_[v] ? 0 : 1);

        int s_mid = max_closure_with_keep_block(m,h,w_scaled,lmid,keep,block,Smid,cuts,cfg_.verbose);

        // 若 Smid 与端点之一相同，尝试往内微调，最多尝试 32 次
        int tries = 0;
        while (tries < 32 && (Smid == Shi_ || Smid == Slo_)) {
            // 选更靠近下端的 λ，以求“扩大” S
            lmid = (Smid == Shi_ ? (lmid + llo)/2 : (lhi + lmid)/2);
            s_mid = max_closure_with_keep_block(m,h,w_scaled,lmid,keep,block,Smid,cuts,cfg_.verbose);
            ++tries;
        }
        if (Smid == Shi_ || Smid == Slo_) {
            // —— 长平台：在子图 D 内做拓扑细化（按权优先） —— //
            // 构子图入度
            vector<int> indeg(n,0);
            for (int i=0; i<m; ++i)
            for (int j=0; j<h; ++j) {
                int v = vid(i,j,h);
                if (!Slo_[v] || Shi_[v] || fixed[v]) continue; // 仅在 D 内
                if (i+1<m) { int u=vid(i+1,j,h); if (Slo_[u] && !Shi_[u] && !fixed[u]) ++indeg[u]; }
                if (j+1<h) { int u=vid(i,j+1,h); if (Slo_[u] && !Shi_[u] && !fixed[u]) ++indeg[u]; }
            }
            struct Key { int w; int i; int j; int v; };
            struct Cmp {
                bool operator()(const Key& A, const Key& B) const {
                    if (A.w != B.w) return A.w < B.w; // 大权优先
                    if (A.i != B.i) return A.i < B.i; // 更靠上优先（i大者先）=> 用 <
                    return A.j < B.j;                 // 更靠右优先
                }
            };
            std::priority_queue<Key, std::vector<Key>, Cmp> pq;
            for (int i=0; i<m; ++i)
            for (int j=0; j<h; ++j) {
                int v = vid(i,j,h);
                if (Slo_[v] && !Shi_[v] && !fixed[v] && indeg[v]==0) {
                    int w = (int)(w_scaled[v]/SCALE); // 原始权
                    pq.push(Key{w,i,j,v});
                }
            }
            while (!pq.empty()) {
                auto k = pq.top(); pq.pop();
                if (fixed[k.v]) continue;
                rank[k.v] = next_rank++; fixed[k.v]=1;
                int i=k.i, j=k.j;
                if (i+1<m) {
                    int u=vid(i+1,j,h);
                    if (Slo_[u] && !Shi_[u] && !fixed[u]) {
                        if (--indeg[u]==0) {
                            int w = (int)(w_scaled[u]/SCALE);
                            pq.push(Key{w,i+1,j,u});
                        }
                    }
                }
                if (j+1<h) {
                    int u=vid(i,j+1,h);
                    if (Slo_[u] && !Shi_[u] && !fixed[u]) {
                        if (--indeg[u]==0) {
                            int w = (int)(w_scaled[u]/SCALE);
                            pq.push(Key{w,i,j+1,u});
                        }
                    }
                }
            }
            return;
        }

        // 正常切入：Smid 在两端之间，递归左右
        refine(lhi, Shi_, lmid, Smid);
        refine(lmid, Smid, llo, Slo_);
    };

    // 启动递归
    refine(lam_hi, Shi, lam_lo, Slo);

    // 有时 platform 细化会一次性定完一大段 rank，若仍有未赋值（极少），按拓扑补齐
    if (next_rank <= n) {
        vector<int> indeg(n,0);
        for (int i=0;i<m;++i) for (int j=0;j<h;++j){
            int v=vid(i,j,h);
            if (i>0) ++indeg[v];
            if (j>0) ++indeg[v];
        }
        struct K2{int w,i,j,v;};
        struct C2{bool operator()(const K2&A,const K2&B)const{
            if (A.w!=B.w) return A.w<B.w;
            if (A.i!=B.i) return A.i<B.i;
            return A.j<B.j;}};
        std::priority_queue<K2,std::vector<K2>,C2> pq;
        for (int i=0;i<m;++i) for (int j=0;j<h;++j){
            int v=vid(i,j,h);
            if (rank[v]==0 && indeg[v]==0){
                int w=(int)(w_scaled[v]/SCALE);
                pq.push(K2{w,i,j,v});
            }
        }
        while(!pq.empty()){
            auto k=pq.top(); pq.pop();
            if (rank[k.v]!=0) continue;
            rank[k.v]=next_rank++;
            int i=k.i,j=k.j;
            if (i+1<m){ int u=vid(i+1,j,h); if (--indeg[u]==0 && rank[u]==0){
                int w=(int)(w_scaled[u]/SCALE); pq.push(K2{w,i+1,j,u});}}
            if (j+1<h){ int u=vid(i,j+1,h); if (--indeg[u]==0 && rank[u]==0){
                int w=(int)(w_scaled[u]/SCALE); pq.push(K2{w,i,j+1,u});}}
        }
    }

    // 回填 y_order
    vector<vector<int>> Y(m, vector<int>(h, 0));
    for (int i=0, v=0; i<m; ++i)
        for (int j=0; j<h; ++j, ++v)
            Y[i][j] = rank[v];

    PTRResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(Y);
    R.x_of_col.assign(h, 0);

    // 成本与核验
    long long Z = 0;
    for (int i=0;i<m;++i)
        for (int j=0;j<h;++j){
            int bw = base_weight(i,j,m,h);
            Z += (long long)bw * (long long)R.y_order[i][j] * cfg_.dV;
        }
    R.cost = Z;
    R.actual_hpwl = hpwl_edges_sum(R.y_order, R.x_of_col, cfg_.dV);
    R.oc_ok = check_global_OC(R.y_order);

    // 唯一占用
    {
        std::unordered_set<uint64_t> occ; occ.reserve((size_t)n*2);
        bool ok = true;
        for (int i=0;i<m;++i)
            for (int j=0;j<h;++j){
                uint64_t x=0ull, yv=(uint64_t)R.y_order[i][j];
                uint64_t key=(x<<32)^yv;
                if (!occ.insert(key).second) ok=false;
            }
        R.unique_ok = ok;
    }
    R.cuts_run = (int)cuts;

    if (cfg_.verbose) {
        cout << "[λ-refine] cuts=" << R.cuts_run
             << ", cost=" << R.cost
             << ", HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok? "OK":"FAIL")
             << ", Unique=" << (R.unique_ok? "OK":"FAIL") << endl;
    }
    return R;
}
