#include "lightPureMcmf.h"
#include <lemon/smart_graph.h>
#include <lemon/network_simplex.h>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cassert>

using Graph = lemon::SmartDigraph;
using std::vector;
using std::cout;
using std::endl;

// 计算相邻边 HPWL（单列验证：水平代价仅取y差，x固定为0）
long long lightPureMcmf::hpwl_sum(const vector<vector<long long>>& y,
                                  const vector<long long>& /*x*/,
                                  long long dV) {
    int m = (int)y.size(), h = (int)y[0].size();
    long long S = 0;
    // 水平邻边
    for (int i = 0; i < m; ++i)
        for (int j = 0; j + 1 < h; ++j)
            S += std::llabs(y[i][j+1] - y[i][j]) * dV;
    // 竖直邻边
    for (int i = 0; i + 1 < m; ++i)
        for (int j = 0; j < h; ++j)
            S += std::llabs(y[i+1][j] - y[i][j]) * dV;
    return S;
}

// 全局 OC：若 i1<=i2 且 j1<=j2，则 y[i1][j1] <= y[i2][j2]
bool lightPureMcmf::check_OC(const vector<vector<long long>>& y) {
    int m = (int)y.size(), h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
        for (int j1 = 0; j1 < h; ++j1)
            for (int i2 = i1; i2 < m; ++i2)
                for (int j2 = j1; j2 < h; ++j2)
                    if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

LightResult lightPureMcmf::solve(int m, int h) {
    using namespace std;

    // 边权重：四边-1/+1，如你原版定义
    auto weight = [&](int i, int j)->long long {
        long long w = 0;
        if (i == 0)     w -= 1;
        if (i == m - 1) w += 1;
        if (j == 0)     w -= 1;
        if (j == h - 1) w += 1;
        return w;
    };

    // ---------- 构图 ----------
    Graph g;
    Graph::Node s = g.addNode(), t = g.addNode();
    Graph::NodeMap<long long> supply(g);       // 节点供需
    Graph::ArcMap<int> cap(g);                 // 容量
    Graph::ArcMap<long long> cost(g);          // 费用
    Graph::NodeMap<int> node2id(g);            // 复原用：Node -> state id（O(1)）

    supply[s] = +1;
    supply[t] = -1;

    struct State { vector<int> a; Graph::Node node; int level = 0; };

    // 10x10 状态规模 ~ C(20,10)=184,756；预留更大空间，减少重分配
    vector<State> states;
    states.reserve(250000);

    // 用 vector<int> 作为 key（比字符串省内存；比 uint64 编码更稳妥，避免溢出）
    struct VecHash {
        size_t operator()(const vector<int>& a) const noexcept {
            // 64-bit 混合哈希（FNV-like + splitmix64）
            uint64_t h = 1469598103934665603ull;
            for (int v : a) {
                uint64_t x = static_cast<uint64_t>(v) + 0x9e3779b97f4a7c15ull;
                h ^= x + (h << 6) + (h >> 2);
            }
            return static_cast<size_t>(h);
        }
    };
    std::unordered_map<vector<int>, int, VecHash> id_of;
    id_of.reserve(250000);

    auto add_state = [&](const vector<int>& a, int lvl) -> int {
        auto it = id_of.find(a);
        if (it != id_of.end()) return it->second;
        State st{a, g.addNode(), lvl};
        supply[st.node] = 0;
        int id = (int)states.size();
        states.push_back(std::move(st));
        id_of.emplace(states.back().a, id);
        node2id[states.back().node] = id;
        return id;
    };

    const int N = m * h;
    vector<int> a0(m, 0);
    int id0 = add_state(a0, 0);

    vector<int> frontier{ id0 }, next; // 仅保留当前层与下一层（节省内存）
    next.reserve(1 << 15);

    for (int lvl = 0; lvl < N; ++lvl) {
        next.clear();
        for (int sid : frontier) {
            const auto& a = states[sid].a;
            // 尝试扩展每一行
            for (int i = 0; i < m; ++i) {
                if (a[i] >= h) continue;
                if (i > 0 && a[i] + 1 > a[i - 1]) continue; // 保持阶梯形非增
                vector<int> b = a;
                int j_new = a[i]; // 扩张前列号即新放置位置（0-based）
                b[i]++;

                int idb = add_state(b, lvl + 1);

                long long w = weight(i, j_new);
                long long c = cfg_.dV * w * (long long)(lvl + 1);

                auto arc = g.addArc(states[sid].node, states[idb].node);
                cap[arc] = 1;
                cost[arc] = c;

                next.push_back(idb);
            }
        }
        // 去重本层产生的 next
        std::sort(next.begin(), next.end());
        next.erase(std::unique(next.begin(), next.end()), next.end());
        if (next.empty()) {
            throw std::runtime_error("Empty next frontier at level " + std::to_string(lvl));
        }
        frontier.swap(next);
    }

    if (frontier.empty()) {
        throw std::runtime_error("Final frontier empty");
    }

    // 源/汇连接
    {
        auto a_arc = g.addArc(s, states[id0].node);
        cap[a_arc] = 1; cost[a_arc] = 0;
    }
    for (int sid : frontier) {
        auto a_arc = g.addArc(states[sid].node, t);
        cap[a_arc] = 1; cost[a_arc] = 0;
    }

    if (cfg_.verbose) {
        cout << "[LightMCMF] states=" << states.size() << " built." << endl;
    }

    // ---------- 最小费用流 ----------
    lemon::NetworkSimplex<Graph, int, long long> ns(g);
    ns.upperMap(cap).costMap(cost).supplyMap(supply);
    if (ns.run() != decltype(ns)::OPTIMAL) {
        throw std::runtime_error("NetworkSimplex failed.");
    }
    long long total_cost = ns.totalCost();

    // ---------- 复原路径（y_order） ----------
    vector<vector<long long>> y(m, vector<long long>(h, 0));
    Graph::Node cur = states[id0].node;
    for (int lvl = 0; lvl < N; ++lvl) {
        Graph::Arc chosen = lemon::INVALID;
        for (Graph::OutArcIt a(g, cur); a != lemon::INVALID; ++a) {
            if (ns.flow(a) > 0) { chosen = a; break; }
        }
        if (chosen == lemon::INVALID) {
            throw std::runtime_error("Path reconstruction failed at level " + std::to_string(lvl));
        }
        Graph::Node nxt = g.target(chosen);

        int sid = node2id[cur];
        int tid = node2id[nxt];

        int sel_i = -1, sel_j = -1;
        for (int i = 0; i < m; ++i) {
            if (states[tid].a[i] == states[sid].a[i] + 1) {
                sel_i = i;
                sel_j = states[sid].a[i]; // 扩张前的列号
                break;
            }
        }
        if (sel_i < 0) throw std::runtime_error("Cannot detect selected cell at level " + std::to_string(lvl));

        y[sel_i][sel_j] = lvl + 1; // rank
        cur = nxt;
    }

    // ---------- 汇总校验 ----------
    LightResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(y);
    R.x_of_col.assign(h, 0);
    R.mcmf_cost = total_cost;
    R.actual_hpwl = hpwl_sum(R.y_order, R.x_of_col, cfg_.dV);
    R.oc_ok = check_OC(R.y_order);

    // 唯一占用（同列）
    {
        std::unordered_set<unsigned long long> occ;
        bool ok = true;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < h; ++j) {
                unsigned long long key = static_cast<unsigned long long>(R.y_order[i][j]);
                if (!occ.insert(key).second) ok = false;
            }
        R.unique_ok = ok;
    }

    if (cfg_.verbose) {
        cout << "  [verify] cost=" << R.mcmf_cost
             << ", HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok ? "OK" : "FAIL")
             << ", Unique=" << (R.unique_ok ? "OK" : "FAIL") << endl;
    }
    return R;
}
