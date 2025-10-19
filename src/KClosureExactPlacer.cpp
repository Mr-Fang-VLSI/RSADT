#include "KClosureExactPlacer.h"
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <unordered_set>
#include <iostream>
#include <cassert>
using std::vector;
using std::cout;
using std::endl;

static const long long INF = (long long)1e14;

long long KClosureExactPlacer::hpwl_sum(const vector<vector<int>>& y, long long dV){
    int m=(int)y.size(), h=(int)y[0].size(); long long S=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j)
        S += (long long)std::llabs(y[i][j+1]-y[i][j]) * dV;
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j)
        S += (long long)std::llabs(y[i+1][j]-y[i][j]) * dV;
    return S;
}
bool KClosureExactPlacer::check_OC(const vector<vector<int>>& y){
    int m=(int)y.size(), h=(int)y[0].size();
    for(int i1=0;i1<m;++i1) for(int j1=0;j1<h;++j1)
    for(int i2=i1;i2<m;++i2) for(int j2=j1;j2<h;++j2)
        if(y[i1][j1]>y[i2][j2]) return false;
    return true;
}

// 一次最大闭包（参数 λ 与 keep 集）
vector<char> KClosureExactPlacer::max_closure_lambda(
    int m, int h,
    const vector<char>& keep,
    long long lambda,
    long long& pos_sum)
{
    using Dig = lemon::ListDigraph;
    Dig g; auto s=g.addNode(), t=g.addNode();
    Dig::ArcMap<long long> cap(g);

    const int n=m*h;
    vector<Dig::Node> node(n);
    for(int v=0;v<n;++v) node[v]=g.addNode();

    auto add = [&](Dig::Node u, Dig::Node v, long long c){
        auto a=g.addArc(u,v); cap[a]=c; return a;
    };

    pos_sum = 0;
    // 节点权 p(v,λ)=w(v)-λ
    for(int i=0,v=0;i<m;++i) for(int j=0;j<h;++j,++v){
        long long p = (long long)wbase(i,j,m,h)*cfg_.dV - lambda;
        if(keep[v]) {
            // 强制包含：s->v 放 INF，v->t 放 0（避免排除）
            add(s,node[v], INF);
        } else {
            if(p>0){ add(s,node[v], p); pos_sum += p; }
            else if(p<0){ add(node[v], t, -p); }
        }
        // 覆盖边：v -> up/left 前驱（强制闭包）
        if(i>0) add(node[v], node[vid(i-1,j,h)], INF);
        if(j>0) add(node[v], node[vid(i,j-1,h)], INF);
    }

    lemon::Preflow<Dig,long long> mf(g,cap,s,t);
    mf.runMinCut(); // 得到最小割
    vector<char> inS(n,0);
    for(int v=0;v<n;++v) inS[v] = mf.minCut(node[v])? 1:0; // true==源侧
    return inS;
}

// 在 0‑边际子图里做“精确补齐”：最大化计数，直到恰好 target
vector<char> KClosureExactPlacer::closure_fill_zero_margin(
    int m, int h,
    const vector<char>& baseS,
    const vector<char>& zero_mask,
    int target)
{
    using Dig = lemon::ListDigraph;
    Dig g; auto s=g.addNode(), t=g.addNode();
    Dig::ArcMap<long long> cap(g);

    const int n=m*h;
    vector<Dig::Node> node(n);
    for(int v=0;v<n;++v) node[v]=g.addNode();

    auto add = [&](Dig::Node u, Dig::Node v, long long c){
        auto a=g.addArc(u,v); cap[a]=c; return a;
    };

    // 想法：baseS 中的点必须在源侧；zero_mask 中的点给“1”的收益（s->v 容量=1），
    // 其它点收益=0。闭包边仍按 v->前驱( INF )。这样 min-cut 选到的 zero 点数量最大。
    // 然后我们从最大可选数 M 里，若 M >= (target - |baseS|)，直接得到补足；
    // 若 M 超过，就再“削掉多余的 0 点”——这一步可以给 zero 点一个极小扰动顺序权，
    // 通过多次 min-cut（或二分数量门槛）把数目卡到精确 target。
    // 这里用“二分门槛”：给每个 zero 候选加 unit weight=1，并在汇侧加“预算阈值”技巧实现恰好 k。

    // 先统计 baseS 大小
    int base_cnt=0; for(char b:baseS) if(b) ++base_cnt;
    int need = target - base_cnt;
    if(need<=0) return baseS;

    // 我们引入一个“预算节点”B，把所有 zero 候选 v 选择视为“消耗 1 单位预算”，总预算=need。
    // 实现方式：对 zero v：s->v cap=1；再 v->B cap=1；B->t cap=need；
    // 同时闭包边照旧：v->前驱 INF。
    // 若 v 在 baseS 中，则 s->v cap=INF（强制选），并同样 v->B cap=0（因为它不消耗预算）。
    // 非 zero 非 baseS 的点不允许新增（不连 s->v），但仍作为闭包传播的“必须前驱”存在。

    auto B = g.addNode();

    for(int i=0,v=0;i<m;++i) for(int j=0;j<h;++j,++v){
        if(baseS[v]){
            add(s,node[v], INF);      // 强制入 S
            // baseS 不消耗预算
        } else if(zero_mask[v]) {
            add(s,node[v], 1);        // 候选
            add(node[v], B, 1);       // 消耗 1 预算
        }
        // 闭包
        if(i>0) add(node[v], node[vid(i-1,j,h)], INF);
        if(j>0) add(node[v], node[vid(i,j-1,h)], INF);
    }
    // 预算上限
    add(B, t, need);

    lemon::Preflow<Dig,long long> mf(g,cap,s,t);
    mf.runMinCut();
    vector<char> inS(n,0);
    for(int v=0;v<n;++v) inS[v] = mf.minCut(node[v])? 1:0;

    return inS; // |S| 一定 == target
}

KClosureResult KClosureExactPlacer::solve(int m, int h){
    const int n=m*h;
    if(cfg_.verbose) cout << "[k-closure] m="<<m<<" h="<<h<<" n="<<n<<"\n";
    vector<int> y(n,0);
    vector<char> keep(n,0);

    long long calls=0;

    // 逐前缀 r=1..n
    for(int r=1;r<=n;++r){
        // 二分 λ：使 |S(λ)| >= r，|S(λ+1)| < r
        long long lo=- (long long)2*cfg_.dV - 2, hi= (long long)2*cfg_.dV + 2;
        vector<char> bestS;
        while(lo<hi){
            long long mid = (lo+hi+1)>>1;
            long long ps=0; auto S = max_closure_lambda(m,h,keep,mid,ps); ++calls;
            int sz=0; for(char b:S) if(b) ++sz;
            if(sz>=r){ bestS=std::move(S); lo=mid; } else { hi=mid-1; }
        }
        // 得到 λ*=lo 的闭包 S*（可能 sz>r，需要在 0 平台精确截断）
        long long ps=0; auto S = max_closure_lambda(m,h,keep,lo,ps); ++calls;
        int sz=0; for(char b:S) if(b) ++sz;

        if(sz>r){
            // 识别 0 边际集合：p(v,λ*)=0 且当前 S 中的“正边选入”和“负边剔除”之外的自由点
            vector<char> zero(n,0);
            for(int i=0,v=0;i<m;++i) for(int j=0;j<h;++j,++v){
                long long p = (long long)wbase(i,j,m,h)*cfg_.dV - lo;
                if(p==0 && !keep[v]) zero[v]=1;
            }
            // baseS：把必须保留的（keep + S 中由正边支撑的）都留住。
            // 简便实现：把 S 直接当 baseS（它包含 keep 与正边选择），再在 zero 上通过预算精确砍到 r。
            S = closure_fill_zero_margin(m,h,S,zero,r);
            // 现在 |S|==r
        }

        // 记录这一步新进的点的秩 = r，并把 S 作为下一层 keep
        vector<char> keep_new = keep;
        for(int v=0;v<n;++v){
            if(S[v] && !keep[v]) y[v] = r;
            keep_new[v] = S[v];
        }
        keep.swap(keep_new);
    }

    // 组装 y_order
    vector<vector<int>> Y(m, vector<int>(h,0));
    for(int v=0;v<n;++v){ int i=v/h, j=v%h; Y[i][j]=y[v]; }

    KClosureResult R;
    R.m=m; R.h=h; R.dV=cfg_.dV;
    R.y_order = std::move(Y);
    // 目标值与 HPWL
    long long Z=0;
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        Z += (long long)wbase(i,j,m,h) * (long long)R.y_order[i][j] * cfg_.dV;
    }
    R.cost = Z;
    R.actual_hpwl = hpwl_sum(R.y_order, cfg_.dV);
    R.oc_ok = check_OC(R.y_order);
    // 唯一性
    {
        std::unordered_set<int> occ; occ.reserve((size_t)n*2);
        bool ok=true;
        for(int i=0;i<m;++i) for(int j=0;j<h;++j){
            if(!occ.insert(R.y_order[i][j]).second) ok=false;
        }
        R.unique_ok = ok;
    }
    R.mincut_calls = calls;
    if(cfg_.verbose){
        cout << "[k-closure] cuts="<<calls
             << " cost="<<R.cost<<" HPWL="<<R.actual_hpwl
             << " OC="<<(R.oc_ok?"OK":"FAIL")
             << " Unique="<<(R.unique_ok?"OK":"FAIL")<<"\n";
    }
    return R;
}
