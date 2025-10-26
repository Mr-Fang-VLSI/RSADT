#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <queue>
#include <string>
#include <tuple>
#include <vector>
using namespace std;

/*
  Route A: 单次 s-t 最小割 == 最大权闭包
  - 图结构：Edmonds–Karp（纯迭代），避免任何递归。
  - 支持: --verbose, --check-oc, --verify, --progress（progress仅占位）
  - 输入: oc_closure m h weighted
      m, h 定义一个 m x h 的网格（总节点 N=m*h）。
      我这里演示了一种典型“正向（右/下）闭包”的约束：选 (x,y) 则必须同时选 (x+1,y) 和 (x,y+1)（若存在）。
      你可以根据你工程内的实际闭包边生成方式，替换 build_closure_edges(...) 即可。
      weighted=0: 所有节点权重=+1
      weighted=1: 示例里用一个既有正也有负的简单确定性权重场（方便验证）；你可以替换为你的 phi/w 计算。
  - 验证:
      --check-oc=1: 返回解是否满足闭包边约束
      --verify=1: 若 N<=22，暴力枚举验证最优性（安全阈值可改）；N>22 会自动跳过并提示
  - 调试:
      --verbose=1: 输出网络统计、INF、可达集等中间信息；--progress=1: 预留（这里不额外输出）
*/

// ========== 简单命令行解析 ==========
static int getFlagInt(int argc, char** argv, const string& key, int default_val) {
    string k = "--" + key + "=";
    for (int i=1;i<argc;i++) {
        string s = argv[i];
        if (s.rfind(k, 0) == 0) {
            return atoi(s.substr(k.size()).c_str());
        }
    }
    return default_val;
}

// ========== 最大流（Edmonds-Karp，完全非递归） ==========
struct Edge {
    int to;
    int rev;            // 反向边在邻接表中的下标
    long long cap;      // 残量容量
    Edge(int _to,int _rev,long long _cap):to(_to),rev(_rev),cap(_cap){}
};

struct MaxFlow {
    int n;
    vector<vector<Edge>> g;
    MaxFlow(int n=0){init(n);}
    void init(int n_) { n=n_; g.assign(n, {}); }

    void addEdge(int u, int v, long long c) {
        Edge a(v, (int)g[v].size(), c);
        Edge b(u, (int)g[u].size(), 0);
        g[u].push_back(a);
        g[v].push_back(b);
    }

    // Edmonds-Karp: BFS 找增广路
    long long maxflow(int s, int t) {
        long long flow = 0;
        vector<int> pv(n), pe(n); // 记录父节点与通过的边号
        while (true) {
            fill(pv.begin(), pv.end(), -1);
            queue<int> q;
            q.push(s);
            pv[s] = s;
            while (!q.empty() && pv[t]==-1) {
                int u = q.front(); q.pop();
                for (int i=0;i<(int)g[u].size();i++) {
                    Edge &e = g[u][i];
                    if (pv[e.to]==-1 && e.cap>0) {
                        pv[e.to] = u; pe[e.to] = i;
                        q.push(e.to);
                        if (e.to==t) break;
                    }
                }
            }
            if (pv[t]==-1) break; // 无增广
            long long aug = LLONG_MAX;
            for (int v=t; v!=s; v=pv[v]) {
                Edge &e = g[pv[v]][pe[v]];
                aug = min(aug, e.cap);
            }
            for (int v=t; v!=s; v=pv[v]) {
                Edge &e = g[pv[v]][pe[v]];
                Edge &r = g[e.to][e.rev];
                e.cap -= aug; r.cap += aug;
            }
            flow += aug;
        }
        return flow;
    }

    // 残量网络从 s 可达性
    vector<char> reachable_from(int s) const {
        vector<char> vis(n, 0);
        queue<int> q; q.push(s); vis[s]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            for (auto &e: g[u]) {
                if (!vis[e.to] && e.cap>0) {
                    vis[e.to]=1; q.push(e.to);
                }
            }
        }
        return vis;
    }
};

// ========== 将 MWC 转为一次最小割 ==========
struct MWCResult {
    vector<int> S;         // 选中的节点集合
    long long obj = 0;     // 目标函数值 sum(w_i) over S
    long long flow = 0;    // 最小割容量（=最大流）
};

MWCResult max_weight_closure(
    int N,
    const vector<pair<int,int>>& must_edges, // i=>j 约束：选 i 必须选 j
    const vector<long long>& w,
    bool verbose
){
    long long sumPos=0, sumNeg=0;
    for (int i=0;i<N;i++){
        if (w[i]>0) sumPos += w[i];
        else        sumNeg += -w[i];
    }
    long long sumAbs = sumPos + sumNeg;
    long long INF = sumAbs + 1; // 强制闭包的足够大容量

    if (verbose) {
        cout << "[closure] |mask|=" << N
             << " sumPos=" << sumPos
             << " sumNeg=" << sumNeg
             << " sumAbs=" << sumAbs
             << " INF="    << INF << "\n";
    }

    int s = N, t = N+1;
    MaxFlow mf(N+2);
    // 节点权重边
    for (int i=0;i<N;i++){
        if (w[i] > 0) mf.addEdge(s, i, w[i]);
        else if (w[i] < 0) mf.addEdge(i, t, -w[i]);
    }
    // 闭包约束边 i=>j
    for (auto &e: must_edges){
        int i=e.first, j=e.second;
        if (i<0||i>=N||j<0||j>=N) continue;
        mf.addEdge(i, j, INF);
    }

    long long flow = mf.maxflow(s, t);
    vector<char> reach = mf.reachable_from(s);

    MWCResult res;
    res.flow = flow;
    for (int i=0;i<N;i++){
        if (reach[i]) {
            res.S.push_back(i);
            res.obj += w[i];
        }
    }

    if (verbose) {
        cout << "=> |S|=" << res.S.size()
             << " sumW(S)=" << res.obj
             << "  (flow=" << res.flow << ")\n";
    }
    return res;
}

// ========== 构建网格上的“正向（右/下）闭包”约束（示例） ==========
// 选 (x,y) 则必须选 (x+1,y) 与 (x,y+1)（若存在）。这是很常见的一类二维正向/上闭包。
static inline int idx(int x, int y, int m, int h){ return y*m + x; }

vector<pair<int,int>> build_closure_edges_grid_right_down(int m, int h) {
    vector<pair<int,int>> E;
    for (int y=0;y<h;y++){
        for (int x=0;x<m;x++){
            int u = idx(x,y,m,h);
            if (x+1<m) E.emplace_back(u, idx(x+1,y,m,h)); // 右
            if (y+1<h) E.emplace_back(u, idx(x,y+1,m,h)); // 下
        }
    }
    return E;
}

// ========== 构建权重（示例） ==========
// 根据 weighted 标志演示两种：
//  - weighted=0: 所有节点 +1
//  - weighted=1: 一个确定性的混合正负权（不随机，便于复现实验/验证）
//   你可以把这里替换成你原来的 phi/w 计算逻辑。
vector<long long> build_weights_example(int m, int h, int weighted) {
    int N = m*h;
    vector<long long> w(N, 1);
    if (weighted==0) {
        // 纯 +1
        fill(w.begin(), w.end(), 1);
    } else {
        // 既有正也有负：例如 w(x,y)=3-(x%4) + ((y%3==0)?1:0) 再把一部分设负
        for (int y=0;y<h;y++){
            for (int x=0;x<m;x++){
                long long val = 3 - (x % 4) + ((y % 3)==0 ? 1 : 0); // 3,2,1,0 循环 + 少量偏置
                if ((x+y)%5==0) val -= 4;  // 每 5 个打一个负冲击
                w[idx(x,y,m,h)] = val;
            }
        }
    }
    return w;
}

// ========== 校验解是否满足闭包 ==========
bool check_closure_feasible(int N, const vector<pair<int,int>>& must_edges, const vector<char>& inS, bool verbose) {
    for (auto &e: must_edges){
        int i=e.first, j=e.second;
        if (inS[i] && !inS[j]) {
            if (verbose) {
                cerr << "[check-oc] violated: choose " << i << " requires choose " << j << "\n";
            }
            return false;
        }
    }
    return true;
}

// ========== 小规模暴力验证（可选） ==========
struct BruteRes { long long best = LLONG_MIN; vector<int> S; };

BruteRes brute_force_mwc(int N, const vector<pair<int,int>>& must_edges, const vector<long long>& w, bool verbose) {
    BruteRes br;
    vector<int> req_out(N, 0);
    vector<vector<int>> succ(N);
    for (auto &e: must_edges) succ[e.first].push_back(e.second);

    auto closed_ok = [&](uint64_t mask)->bool{
        for (int i=0;i<N;i++){
            if (!(mask>>i & 1ULL)) continue;
            for (int j: succ[i]) {
                if (!(mask>>j & 1ULL)) return false;
            }
        }
        return true;
    };

    uint64_t total = (N>=63)? 0ULL : (1ULL<<N);
    for (uint64_t s=0; s<total; ++s){
        if (!closed_ok(s)) continue;
        long long val=0;
        for (int i=0;i<N;i++) if (s>>i & 1ULL) val+=w[i];
        if (val>br.best){ br.best=val; br.S.clear(); for(int i=0;i<N;i++) if (s>>i & 1ULL) br.S.push_back(i); }
    }
    if (verbose) {
        if (total==0) cerr << "[verify] N too large for brute force\n";
        else cerr << "[verify] brute enumerated " << total << " subsets\n";
    }
    return br;
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " m h weighted"
             << " [--verify=0|1] [--check-oc=0|1] [--verbose=0|1] [--progress=0|1]\n";
        return 1;
    }
    int m = atoi(argv[1]);
    int h = atoi(argv[2]);
    int weighted = atoi(argv[3]);
    int verify   = getFlagInt(argc, argv, "verify",   0);
    int check_oc = getFlagInt(argc, argv, "check-oc", 0);
    int verbose  = getFlagInt(argc, argv, "verbose",  0);
    int progress = getFlagInt(argc, argv, "progress", 0); // 这里不额外使用，仅占位保持接口

    int N = m*h;
    cout << "[Closure] m=" << m << " h=" << h << " N=" << N
         << " (closure exact, weighted=" << weighted << ")\n";

    // 1) 构建闭包边（你可以替换成本项目真实的闭包边）
    vector<pair<int,int>> must_edges = build_closure_edges_grid_right_down(m,h);

    // 2) 构建权重（你可以替换成你的 phi/w）
    vector<long long> w = build_weights_example(m,h,weighted);

    // 3) 跑一次最大权闭包（单次 s-t 最大流）
    auto res = max_weight_closure(N, must_edges, w, verbose);

    // 输出解要素
    cout << "[Result] |S|=" << res.S.size() << " Obj=" << res.obj << "\n";
    if (verbose) {
        cout << "[Result] S indices: ";
        for (size_t i=0;i<res.S.size();i++){
            cout << res.S[i] << (i+1==res.S.size()?'\n':' ');
        }
    }

    // 4) 可选：检查闭包可行
    if (check_oc) {
        vector<char> inS(N, 0);
        for (int i: res.S) inS[i]=1;
        bool ok = check_closure_feasible(N, must_edges, inS, verbose);
        cout << "[check-oc] " << (ok? "OK" : "FAILED") << "\n";
        if (!ok) return 2;
    }

    // 5) 可选：小规模暴力验证
    if (verify) {
        if (N <= 22) {
            auto br = brute_force_mwc(N, must_edges, w, verbose);
            cout << "[verify] brute best=" << br.best << "  flow_obj=" << res.obj << "\n";
            if (br.best != res.obj) {
                cerr << "[verify] MISMATCH! brute says " << br.best
                     << " but flow gives " << res.obj << "\n";
                return 3;
            } else {
                cout << "[verify] MATCH\n";
            }
        } else {
            cout << "[verify] N=" << N << " too large; skip brute force (set a larger threshold if needed)\n";
        }
    }

    return 0;
}
