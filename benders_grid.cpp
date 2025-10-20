#include <bits/stdc++.h>
using namespace std;

/*
  GridBenders: Benders/边界变量主问题 + 差分约束可行性从问题（oracle）
  问题：
    - m x h 网格，每个格点一个变量 y(i,j) ∈ {1..n}, n = m*h
    - OC（行/列严格递增）: 邻接 (右/下) 边满足 y_v - y_u ≥ 1
    - Max-T: 邻接 (右/下) 边满足 y_v - y_u ≤ T
    - AllDifferent: { y(i,j) } 恰好是 {1..n} 的一个排列（单槽一元）
    - 目标（只依赖边界，见 R-SAD Lemma 1 / 式(2)）：min Σ w_u * y_u ，w∈{-1,0,+1}

  方法概要：
    - 主问题：仅枚举边界变量（O(m+h)个），保持 OC & Max-T & [1..n]，并禁止边界之间重复槽位。
      回溯 + 差分传播求解，利用目标符号优先选取边界变量的极值，配合下界剪枝。
    - 从问题（oracle）：
        1) 把已定的边界值注入全图，做差分约束闭包传播，得到每个内部点 [L,U] 区间。
           传播规则（对每条有向边 u→v）：
              L[v] ≥ L[u] + 1;        U[u] ≤ U[v] - 1;     // OC
              U[v] ≤ U[u] + T;        L[u] ≥ L[v] - T;     // Max-T
           如有 L>U 则不可行；
        2) 若可行，做“无空隙”的一机调度（每时隙恰好安排一个点，时间=层级）：
           · 已定边界当作已安排（占用其 y）。
           · 仅当一个内部点的所有前驱都已安排才“就绪”，其可用窗口为
             [ max(L, max_pred(y_pred+1)), min(U, min_pred(y_pred+T)) ]。
           · 逐时隙 t=1..n 递增，若 t 没被边界占用，则从“就绪且 L≤t≤U”的点中
             选截止期（最小 U）最小者放置；必要时回溯，保证 AllDifferent。
        3) 成功则得到全局 y 的整数延拓；目标只看边界即可。
*/

struct GridBenders {
    struct Node {
        int i, j, id;
        vector<int> preds, succs;
        bool boundary = false;
    };

    int m, h, T, n;
    vector<Node> nodes;
    vector<vector<int>> preds, succs;
    vector<int> weights;           // w(i,j) ∈ {-1,0,+1}
    vector<int> boundary_ids;      // 边界变量的节点 id
    vector<int> best_y;            // 最优解
    long long best_obj = (1LL<<62);

    GridBenders(int m_, int h_, int T_) : m(m_), h(h_), T(T_) {
        n = m*h;
        build_graph();
        build_weights();
        collect_boundary();
    }

    // id <-> (i,j)
    inline int id_of(int i, int j) const { return (i-1)*h + (j-1); }
    inline pair<int,int> ij_of(int id) const { return { id / h + 1, id % h + 1 }; }

    void build_graph() {
        nodes.resize(n);
        preds.assign(n, {});
        succs.assign(n, {});
        for (int i=1; i<=m; ++i) {
            for (int j=1; j<=h; ++j) {
                int id = id_of(i,j);
                nodes[id].i = i; nodes[id].j = j; nodes[id].id = id;
                // 边界判定
                if (i==1 || i==m || j==1 || j==h) nodes[id].boundary = true;
                // 邻接（右/下）有向边：u -> v
                if (i+1<=m) { // down
                    int v = id_of(i+1,j);
                    succs[id].push_back(v);
                    preds[v].push_back(id);
                }
                if (j+1<=h) { // right
                    int v = id_of(i,j+1);
                    succs[id].push_back(v);
                    preds[v].push_back(id);
                }
            }
        }
        for (int id=0; id<n; ++id) {
            nodes[id].preds = preds[id];
            nodes[id].succs = succs[id];
        }
    }

    void build_weights() {
        weights.assign(n, 0);
        for (int id=0; id<n; ++id) {
            auto [i,j] = ij_of(id);
            int w = 0;
            if (i==m) w += 1;      // bottom row
            if (j==h) w += 1;      // right col
            if (i==1) w -= 1;      // top row
            if (j==1) w -= 1;      // left col
            weights[id] = w;       // 角点自然合并（例如 (1,1) 得 -2, (m,h) 得 +2）
        }
    }

    void collect_boundary() {
        boundary_ids.clear();
        boundary_ids.reserve(2*(m+h));
        vector<char> seen(n, 0);
        auto push_if = [&](int id){
            if (!seen[id]) { seen[id]=1; boundary_ids.push_back(id); }
        };
        // 统一次序：上→右→下→左（避免角点重复）
        for (int j=1; j<=h; ++j) push_if(id_of(1,j));
        for (int i=2; i<=m; ++i) push_if(id_of(i,h));
        for (int j=h-1; j>=1; --j) push_if(id_of(m,j));
        for (int i=m-1; i>=2; --i) push_if(id_of(i,1));
    }

    // 差分约束闭包传播：给定部分已定 y，求每个点的 [L,U]，并检测 L<=U
    bool propagate_bounds(const vector<int>& y_fixed, vector<int>& L, vector<int>& U) const {
        L.assign(n, 1);
        U.assign(n, n);
        // 赋予固定边界（或其它）值
        for (int id=0; id<n; ++id) if (y_fixed[id] != -1) { L[id] = U[id] = y_fixed[id]; }
        bool changed = true;
        int iter = 0, iter_limit = 4*n; // 足够的上界
        while (changed && iter++ < iter_limit) {
            changed = false;
            // 前向扫描
            for (int id=0; id<n; ++id) {
                for (int v: succs[id]) {
                    // OC: y[v] >= y[id] + 1  ==> L[v] >= L[id]+1 ; U[id] <= U[v]-1
                    int nv = L[id] + 1;
                    if (L[v] < nv) { L[v] = nv; changed = true; }
                    int nu = U[v] - 1;
                    if (U[id] > nu) { U[id] = nu; changed = true; }
                    // Max-T: y[v] <= y[id] + T, y[id] >= y[v] - T
                    int uv = U[id] + T;
                    if (U[v] > uv) { U[v] = uv; changed = true; }
                    int lu = L[v] - T;
                    if (L[id] < lu) { L[id] = lu; changed = true; }
                }
            }
            // 反向扫描
            for (int id=n-1; id>=0; --id) {
                for (int v: succs[id]) {
                    int nv = L[id] + 1;
                    if (L[v] < nv) { L[v] = nv; changed = true; }
                    int nu = U[v] - 1;
                    if (U[id] > nu) { U[id] = nu; changed = true; }
                    int uv = U[id] + T;
                    if (U[v] > uv) { U[v] = uv; changed = true; }
                    int lu = L[v] - T;
                    if (L[id] < lu) { L[id] = lu; changed = true; }
                }
            }
            // 截断到 [1..n]
            for (int id=0; id<n; ++id) {
                if (L[id] < 1) { L[id]=1; }
                if (U[id] > n) { U[id]=n; }
            }
        }
        // 判断可行
        for (int id=0; id<n; ++id) {
            if (L[id] > U[id]) return false;
        }
        return true;
    }

    // 计算目标函数（只看边界）
    long long objective_by_weights(const vector<int>& y) const {
        long long obj = 0;
        for (int id: boundary_ids) {
            if (y[id] == -1) return (1LL<<61); // 未完全赋值不该来这里
            obj += 1LL * weights[id] * y[id];
        }
        return obj;
    }

    // 目标下界（用于主问题剪枝）：已定部分取真值，未定部分 w>=0 取 L，w<0 取 U
    long long objective_lower_bound(const vector<int>& y, const vector<int>& L, const vector<int>& U) const {
        long long lb = 0;
        for (int id: boundary_ids) {
            if (y[id] != -1) lb += 1LL * weights[id] * y[id];
            else {
                if (weights[id] >= 0) lb += 1LL * weights[id] * L[id];
                else                  lb += 1LL * weights[id] * U[id];
            }
        }
        return lb;
    }

    // 内部调度（精确，回溯）：给定边界 y 与 [L,U]，构造 AllDifferent 的全局延拓
    bool schedule_interior(const vector<int>& y_fixed, const vector<int>& L, const vector<int>& U, vector<int>& y_out) const {
        y_out = y_fixed;
        vector<int> indeg(n, 0);
        for (int v=0; v<n; ++v) indeg[v] = (int)preds[v].size();

        // usedTime[t]=1 表示槽位 t 已被边界占用
        vector<char> used(n+1, 0);
        for (int id=0; id<n; ++id) if (y_out[id] != -1) {
            if (used[y_out[id]]) return false; // 边界之间不得重复
            used[y_out[id]] = 1;
        }

        // 动态窗口：随着前驱被安排而收紧
        vector<int> dynL = L, dynU = U;
        vector<int> remPred = indeg;
        vector<char> scheduled(n, 0);

        // 将已经安排（边界）的影响传播到其后继的动态窗口
        for (int u=0; u<n; ++u) if (y_out[u] != -1) {
            scheduled[u] = 1;
            for (int v: succs[u]) {
                remPred[v]--;
                dynL[v] = max(dynL[v], y_out[u] + 1);
                dynU[v] = min(dynU[v], y_out[u] + T);
            }
        }

        // ready 集合：所有 remPred==0 且未安排的点，但我们每个时隙只考虑 dynL<=t<=dynU 的
        vector<int> order; order.reserve(n);
        function<bool(int)> dfs = [&](int t)->bool {
            // 找到下一个未被占用的时隙
            while (t<=n && used[t]) ++t;
            if (t>n) return true; // 所有槽位均已填满

            // 收集此时可放置的“就绪”内部点
            vector<int> cand;
            cand.reserve(8);
            for (int v=0; v<n; ++v) {
                if (scheduled[v]) continue;
                if (remPred[v]!=0) continue;
                // 只允许内部点在这里放置（边界已安排）
                if (nodes[v].boundary) continue;
                if (dynL[v] <= t && t <= dynU[v]) cand.push_back(v);
            }
            if (cand.empty()) return false; // 此时隙无法安放任何点 -> 不可行

            // 最小截止期（dynU）优先，减少回溯
            sort(cand.begin(), cand.end(), [&](int a, int b){
                if (dynU[a] != dynU[b]) return dynU[a] < dynU[b];
                return a < b;
            });

            for (int v: cand) {
                // 安排 v @ t
                scheduled[v] = 1;
                y_out[v] = t;
                used[t] = 1;
                // 备份后继的状态以便回溯
                vector<pair<int,pair<int,int>>> backup; backup.reserve(succs[v].size());
                for (int w: succs[v]) {
                    backup.push_back({w, {dynL[w], dynU[w]}});
                    remPred[w]--;
                    dynL[w] = max(dynL[w], t+1);
                    dynU[w] = min(dynU[w], t+T);
                }
                // 就绪时校验窗口不空
                bool ok = true;
                for (auto &bk: backup) {
                    int w = bk.first;
                    if (remPred[w]==0) {
                        // 到达就绪，进一步与静态 L/U 交汇（虽然 dynL/U 本就从静态起步）
                        dynL[w] = max(dynL[w], L[w]);
                        dynU[w] = min(dynU[w], U[w]);
                        if (dynL[w] > dynU[w]) { ok=false; break; }
                    }
                }
                if (ok && dfs(t+1)) return true;

                // 回溯
                for (auto &bk: backup) {
                    int w = bk.first;
                    remPred[w]++;
                    dynL[w] = bk.second.first;
                    dynU[w] = bk.second.second;
                }
                used[t] = 0;
                y_out[v] = -1;
                scheduled[v] = 0;
            }
            return false;
        };

        // 调度
        bool feasible = dfs(1);
        if (!feasible) return false;

        // 最终校验（健壮性）
        vector<char> seen(n+1, 0);
        for (int id=0; id<n; ++id) {
            if (y_out[id] < 1 || y_out[id] > n) return false;
            if (seen[y_out[id]]) return false;
            seen[y_out[id]] = 1;
        }
        for (int u=0; u<n; ++u) for (int v: succs[u]) {
            int d = y_out[v] - y_out[u];
            if (d < 1 || d > T) return false;
        }
        return true;
    }

    // 选择下一个边界变量（最小可行域优先，次序以 |w| 降序）
    int pick_next_boundary(const vector<int>& y, const vector<int>& L, const vector<int>& U, const vector<char>& used) const {
        int best = -1;
        long long bestKey = (1LL<<60);
        for (int id: boundary_ids) {
            if (y[id] != -1) continue;
            // 粗估可行域大小（剔除已用槽位的数量）
            int lo = L[id], hi = U[id];
            int width = max(0, hi - lo + 1);
            // 更“紧”的优先，其次 |w| 越大越优（更影响目标）
            long long key = 1LL*width*1000 - 1LL*abs(weights[id]);
            if (key < bestKey) { bestKey = key; best = id; }
        }
        return best;
    }

    // 主问题回溯 + 剪枝 + 从问题 oracle
    void dfs_master(vector<int>& y, vector<char>& used) {
        // 差分传播
        vector<int> L, U;
        if (!propagate_bounds(y, L, U)) return;

        // 目标下界剪枝
        long long lbObj = objective_lower_bound(y, L, U);
        if (lbObj >= best_obj) return;

        // 是否已完成边界赋值
        bool all_assigned = true;
        for (int id: boundary_ids) if (y[id] == -1) { all_assigned = false; break; }

        if (all_assigned) {
            // 边界不可重复
            for (int id: boundary_ids) if (y[id] != -1) {
                if (used[y[id]] == 0) { /*shouldn't happen*/ }
            }
            // 从问题：尝试构造 AllDifferent 的全局整数延拓
            vector<int> y_full = y;
            if (!schedule_interior(y, L, U, y_full)) return;
            long long obj = objective_by_weights(y_full);
            if (obj < best_obj) { best_obj = obj; best_y = y_full; }
            return;
        }

        // 挑变量
        int var = pick_next_boundary(y, L, U, used);
        if (var == -1) return;

        int lo = L[var], hi = U[var];
        if (lo > hi) return;

        // 按权重符号决定遍历方向：w>=0 低到高（越低越好），w<0 高到低（越高越好）
        vector<int> dom;
        dom.reserve(hi-lo+1);
        if (weights[var] >= 0) {
            for (int t=lo; t<=hi; ++t) dom.push_back(t);
        } else {
            for (int t=hi; t>=lo; --t) dom.push_back(t);
        }

        // 枚举取值（剔除已用槽位）
        for (int val: dom) {
            if (used[val]) continue; // 边界间严格 AllDifferent
            // 早期一致性：检查与已定相邻边界的 OC/MaxT 在传播里已反映，这里不重复
            y[var] = val; used[val] = 1;
            dfs_master(y, used);
            used[val] = 0; y[var] = -1;
        }
    }

    bool solve() {
        best_obj = (1LL<<62);
        best_y.assign(n, -1);
        vector<int> y(n, -1);
        vector<char> used(n+1, 0);

        // 可加一条全局必要割：y(m,h)-y(1,1) ∈ [ (m-1)+(h-1), ((m-1)+(h-1))*T ]，有利于早剪枝
        // 这里不显式加为等式，依赖传播即可（边界回溯时自适应满足）

        dfs_master(y, used);
        return best_obj < (1LL<<62);
    }

    void print_solution() const {
        if (best_obj >= (1LL<<62)) {
            cout << "No feasible solution found.\n";
            return;
        }
        cout << "Best objective (boundary linear form) = " << best_obj << "\n";
        // 同时计算 HPWL 纵向部分（式(2)）校验
        long long hpwl_v = 0;
        for (int j=1; j<=h; ++j) {
            hpwl_v += best_y[id_of(m,j)] - best_y[id_of(1,j)];
        }
        for (int i=1; i<=m; ++i) {
            hpwl_v += best_y[id_of(i,h)] - best_y[id_of(i,1)];
        }
        cout << "HPWL vertical by Eq.(2) = " << hpwl_v << "\n";

        // 打印 y 矩阵
        cout << "y grid (m="<<m<<", h="<<h<<", T="<<T<<"):\n";
        for (int i=1; i<=m; ++i) {
            for (int j=1; j<=h; ++j) {
                int id = id_of(i,j);
                cout << setw(4) << best_y[id];
            }
            cout << "\n";
        }

        // 快速校验
        bool ok = true;
        vector<char> seen(n+1, 0);
        for (int id=0; id<n; ++id) {
            int yv = best_y[id];
            if (yv<1 || yv>n || seen[yv]) ok=false; else seen[yv]=1;
        }
        for (int u=0; u<n; ++u) for (int v: succs[u]) {
            int d = best_y[v] - best_y[u];
            if (d<1 || d> T) ok=false;
        }
        if (!ok) cout << "WARNING: verification failed (should not happen).\n";
    }
};

// ---- demo main ----
int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int m = 4, h = 4, T = 3; // 默认小规模演示
    if (argc >= 3) {
        m = atoi(argv[1]);
        h = atoi(argv[2]);
    }
    if (argc >= 4) T = atoi(argv[3]);

    GridBenders solver(m, h, T);
    bool ok = solver.solve();
    if (!ok) {
        cout << "No feasible solution for m="<<m<<", h="<<h<<", T="<<T<<"\n";
        return 0;
    }
    solver.print_solution();
    return 0;
}
