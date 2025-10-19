#include "PureMcmfColumnPlacer.h"
#include <lemon/list_graph.h>
#include <lemon/network_simplex.h>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>

using Graph = lemon::ListDigraph;
using std::vector;
using std::string;
using std::cout;
using std::endl;

static string key_of(const vector<int>& a) {
    std::ostringstream oss;
    for (size_t i = 0; i < a.size(); ++i) {
        if (i) oss << ',';
        oss << a[i];
    }
    return oss.str();
}

// HPWL（四邻边）
long long PureMcmfColumnPlacer::hpwl_edges_sum(const vector<vector<long long>>& y,
                                               const vector<long long>& x,
                                               long long dV) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long S = 0;
    // 水平邻边（dx=0，因为全在一根 site 列）
    for (int i = 0; i < m; ++i)
        for (int j = 0; j+1 < h; ++j) {
            long long dy = std::llabs(y[i][j+1] - y[i][j]) * dV;
            S += dy;
        }
    // 竖向邻边
    for (int i = 0; i+1 < m; ++i)
        for (int j = 0; j < h; ++j) {
            long long dy = std::llabs(y[i+1][j] - y[i][j]) * dV;
            S += dy;
        }
    return S;
}

// 全局 OC：若 i1<=i2 且 j1<=j2，则 y[i1][j1] <= y[i2][j2]
bool PureMcmfColumnPlacer::check_global_OC(const vector<vector<long long>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
    for (int j1 = 0; j1 < h; ++j1)
    for (int i2 = i1; i2 < m; ++i2)
    for (int j2 = j1; j2 < h; ++j2)
        if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

PureResult PureMcmfColumnPlacer::solve(int m, int h) {
    // 约定：i=0 是底行，i=m-1 是顶行；j=0 是左列，j=h-1 是右列。
    // 权重 w(i,j) = (top?+1 : bottom?-1 : 0) + (right?+1 : left?-1 : 0)
    auto weight = [&](int i, int j)->long long {
        long long w = 0;
        if (i == 0)     w -= 1;
        if (i == m-1)   w += 1;
        if (j == 0)     w -= 1;
        if (j == h-1)   w += 1;
        return w;
    };

    // ---------- 构建“阶梯前沿”状态图 ----------
    // 状态 a[0..m-1]：第 i 行（自底向上）的“已选列数”，满足 0<=a[m-1]<=...<=a[0]<=h。
    // 从状态 a 选择一行 i，把该行从 a[i] 增到 a[i]+1 合法当且仅当：
    //  1) a[i] < h；2) (i==0) 或 (a[i] + 1 <= a[i-1])。新加入的单元是 (i, j=a[i]+1)（1-based 为 a[i]+1）。
    const int N = m * h;

    Graph g;
    Graph::Node s = g.addNode();
    Graph::Node t = g.addNode();
    lemon::ListDigraph::NodeMap<long long> supply(g);
    lemon::ListDigraph::ArcMap<int>        cap(g);
    lemon::ListDigraph::ArcMap<long long>  cost(g);

    supply[s] = +1;
    supply[t] = -1;

    // 各“状态”节点
    struct StateInfo { vector<int> a; int level = 0; Graph::Node node; };
    vector<StateInfo> states; states.reserve(20000);
    std::unordered_map<string,int> id_of; id_of.reserve(20000);

    auto add_state = [&](const vector<int>& a, int level)->int {
        string key = key_of(a);
        auto it = id_of.find(key);
        if (it != id_of.end()) return it->second;
        StateInfo info; info.a = a; info.level = level; info.node = g.addNode();
        supply[info.node] = 0;
        int id = (int)states.size();
        states.push_back(std::move(info));
        id_of.emplace(std::move(key), id);
        return id;
    };

    // 层序生成（BFS）：level = 已选单元数 = sum(a)
    vector<vector<int>> level_ids(N+1);
    vector<int> a0(m, 0);
    int id0 = add_state(a0, 0);
    level_ids[0].push_back(id0);

    for (int lvl = 0; lvl < N; ++lvl) {
        for (int sid : level_ids[lvl]) {
            const auto& a = states[sid].a;
            // 枚举可扩展行
            for (int i = 0; i < m; ++i) {
                if (a[i] >= h) continue;
                if (i > 0 && a[i] + 1 > a[i-1]) continue;
                vector<int> b = a;
                int j_new0 = a[i];        // 0-based j（扩展前值）
                b[i] = a[i] + 1;
                int idb = add_state(b, lvl+1);
                if (states[idb].level == 0) states[idb].level = lvl+1;
                level_ids[lvl+1].push_back(idb);

                // 新增单元 (i, j=j_new0) 其 rank = lvl+1（1-based）
                long long w = weight(i, j_new0);
                long long c = cfg_.dV * w * (long long)(lvl + 1);

                auto a_arc = g.addArc(states[sid].node, states[idb].node);
                cap[a_arc]  = 1;
                cost[a_arc] = c;
            }
        }
        // 去重（同层可能重复 push 相同状态）
        std::sort(level_ids[lvl+1].begin(), level_ids[lvl+1].end());
        level_ids[lvl+1].erase(std::unique(level_ids[lvl+1].begin(), level_ids[lvl+1].end()), level_ids[lvl+1].end());
    }

    // 源/汇连接
    {
        auto a_arc = g.addArc(s, states[id0].node);
        cap[a_arc] = 1; cost[a_arc] = 0;
    }
    for (int sid : level_ids[N]) {
        auto a_arc = g.addArc(states[sid].node, t);
        cap[a_arc] = 1; cost[a_arc] = 0;
    }

    if (cfg_.verbose) {
        cout << "[PureMCMF] build states=" << states.size()
             << ", arcs=";
        long long A = 0;
        for (int lv = 0; lv < N; ++lv) {
            // 每层出边数近似是可扩展行数的总和
            A += (long long)level_ids[lv].size(); // 粗略
        }
        cout << "≈" << A << " (rough)\n";
    }

    // ---------- 最小费用流 ----------
    lemon::NetworkSimplex<Graph,int,long long> ns(g);
    ns.upperMap(cap).costMap(cost).supplyMap(supply);
    auto status = ns.run();
    if (status != decltype(ns)::OPTIMAL) {
        throw std::runtime_error("NetworkSimplex failed on pure MCMF column placement.");
    }
    long long Z = ns.totalCost();

    // ---------- 复原路径（得到 y_order 1..N） ----------
    vector<vector<long long>> y(m, vector<long long>(h, 0));
    // 先反查每条“状态边”对应的新单元 (i,j)
    // 为节省内存，我们动态二次遍历：逐层从当前状态找 flow=1 的出弧，并由 a→b 的差分确定 (i,j)
    Graph::Node cur = states[id0].node;
    for (int lvl = 0; lvl < N; ++lvl) {
        Graph::Arc chosen = lemon::INVALID;
        for (Graph::OutArcIt a(g, cur); a != lemon::INVALID; ++a) {
            if (ns.flow(a) > 0) { chosen = a; break; }
        }
        if (chosen == lemon::INVALID) {
            throw std::runtime_error("Flow path reconstruction failed.");
        }
        Graph::Node nxt = g.target(chosen);

        // 找出 a->b 哪一行增加了 1（即新放置的 (i, j)）
        int sid = -1, tid = -1;
        // 反查节点 id
        for (size_t k = 0; k < states.size(); ++k) {
            if (states[k].node == cur)  sid = (int)k;
            if (states[k].node == nxt)  tid = (int)k;
        }
        if (sid < 0 || tid < 0) throw std::runtime_error("internal id lookup failed");

        int sel_i = -1, sel_j0 = -1;
        for (int i = 0; i < m; ++i) {
            if (states[tid].a[i] == states[sid].a[i] + 1) {
                sel_i = i;
                sel_j0 = states[sid].a[i]; // 0-based j
                break;
            }
        }
        if (sel_i < 0) throw std::runtime_error("cannot detect selected cell");

        y[sel_i][sel_j0] = lvl + 1; // rank

        cur = nxt;
    }

    // ---------- 汇总与校验 ----------
    PureResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(y);
    R.x_of_col.assign(h, 0);
    R.mcmf_cost = Z;
    R.actual_hpwl = hpwl_edges_sum(R.y_order, R.x_of_col, cfg_.dV);
    R.oc_ok = check_global_OC(R.y_order);

    // 唯一占用 (x,y)
    {
        std::unordered_set<unsigned long long> occ;
        bool ok = true;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < h; ++j) {
                unsigned long long x = 0ULL; // 同一列
                unsigned long long yv = (unsigned long long)(R.y_order[i][j] & 0xffffffffULL);
                unsigned long long key = (x << 32) ^ yv;
                if (!occ.insert(key).second) ok = false;
            }
        R.unique_ok = ok;
    }

    if (cfg_.verbose) {
        cout << "  [verify] MCMF cost=" << R.mcmf_cost
             << ", actual HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok ? "OK" : "FAIL")
             << ", Unique=" << (R.unique_ok ? "OK" : "FAIL") << endl;
    }
    return R;
}
