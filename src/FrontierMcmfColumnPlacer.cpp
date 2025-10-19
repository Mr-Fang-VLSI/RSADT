#include "FrontierMcmfColumnPlacer.h"
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <iostream>
#include <algorithm>  // for std::reverse

using std::vector;
using std::pair;
using std::cout;
using std::endl;

// === 真值 HPWL：四邻边 ===
long long FrontierMcmfColumnPlacer::hpwl_edges_sum(const vector<vector<int>>& y,
                                                   const vector<long long>& x,
                                                   long long dV) {
    (void)x; // 单列场景未用到 x，消除 unused-parameter 警告
    int m = (int)y.size();
    int h = (int)y[0].size();
    long long S = 0;
    // 水平邻边（dx=0）
    for (int i = 0; i < m; ++i)
        for (int j = 0; j + 1 < h; ++j) {
            long long dy = std::llabs((long long)y[i][j+1] - y[i][j]) * dV;
            S += dy;
        }
    // 竖向邻边
    for (int i = 0; i + 1 < m; ++i)
        for (int j = 0; j < h; ++j) {
            long long dy = std::llabs((long long)y[i+1][j] - y[i][j]) * dV;
            S += dy;
        }
    return S;
}

// === 全局 OC：若 i1<=i2 且 j1<=j2，则 y[i1][j1] <= y[i2][j2] ===
bool FrontierMcmfColumnPlacer::check_global_OC(const vector<vector<int>>& y) {
    int m = (int)y.size();
    int h = (int)y[0].size();
    for (int i1 = 0; i1 < m; ++i1)
    for (int j1 = 0; j1 < h; ++j1)
    for (int i2 = i1; i2 < m; ++i2)
    for (int j2 = j1; j2 < h; ++j2)
        if (y[i1][j1] > y[i2][j2]) return false;
    return true;
}

// === 状态压缩：前沿向量 a[0..m-1]，每行已取的列数（0..h），且 a[0]>=a[1]>=...>=a[m-1] ===
// 用 8bit 存每个 a[i]（m<=30,h<=30），打包成 string 作为哈希键，节省内存与拷贝。
struct PackedFrontier {
    std::string bytes; // 长度 m，每个字节= a[i]
    bool operator==(const PackedFrontier& o) const { return bytes == o.bytes; }
};

struct PackedHash {
    size_t operator()(const PackedFrontier& k) const noexcept {
        // 64-bit FNV-1a
        uint64_t h = 1469598103934665603ull;
        for (unsigned char c : k.bytes) {
            h ^= (uint64_t)c;
            h *= 1099511628211ull;
        }
        return (size_t)h;
    }
};

FrontierResult FrontierMcmfColumnPlacer::solve(int m, int h) {
    const int N = m * h;

    // ---------- Dijkstra on implicit DAG (exact shortest path = 1-unit MCMF) ----------
    struct NodeRec {
        long long dist;
        int level;     // 已放置单元数
        int parent_id; // 回溯
        int put_i;     // 本步选择的行
        int put_j;     // 本步选择的列
    };

    // 优先队列项：加入自增序号，避免比较 PackedFrontier
    struct QItem {
        long long dist;
        uint64_t  seq;
        PackedFrontier key;
    };
    struct QCmp {
        bool operator()(const QItem& a, const QItem& b) const {
            if (a.dist != b.dist) return a.dist > b.dist; // 最小堆
            return a.seq  > b.seq;                        // 二级关键字
        }
    };

    // 初始化前沿 a=0
    PackedFrontier start;
    start.bytes.assign(m, (char)0);

    // Dijkstra 容器
    std::priority_queue<QItem, std::vector<QItem>, QCmp> pq;
    uint64_t seqgen = 0;

    std::unordered_map<PackedFrontier, NodeRec, PackedHash> best; best.reserve(1<<20);
    std::unordered_map<PackedFrontier, int, PackedHash> idmap; idmap.reserve(1<<20);
    vector<PackedFrontier> id2key; id2key.reserve(1<<20);

    auto push_state = [&](const PackedFrontier& k,
                          long long ndist, int level,
                          int parent_id, int put_i, int put_j) {
        auto it = best.find(k);
        if (it == best.end() || ndist < it->second.dist) {
            best[k] = NodeRec{ndist, level, parent_id, put_i, put_j};
            pq.push(QItem{ndist, seqgen++, k});
            if (idmap.find(k) == idmap.end()) {
                int nid = (int)id2key.size();
                idmap.emplace(k, nid);
                id2key.push_back(k);
            }
        }
    };

    push_state(start, 0, 0, -1, -1, -1);

    PackedFrontier goal;
    long long Z = -1;

    size_t expanded = 0;

    while (!pq.empty()) {
        auto item = pq.top(); pq.pop();
        long long cdist = item.dist;
        PackedFrontier key = item.key;

        auto it = best.find(key);
        if (it == best.end() || cdist != it->second.dist) continue; // 过期
        const NodeRec cur = it->second;
        const int level = cur.level;
        if (level == N) {
            goal = key;
            Z = cdist;
            break;
        }

        // 生成后继：选择一行 i，使 a[i]<h 且 (i==0 || a[i]+1 <= a[i-1])
        const std::string& a = key.bytes;
        for (int i = 0; i < m; ++i) {
            unsigned char ai = (unsigned char)a[i];
            if (ai >= (unsigned char)h) continue;
            if (i > 0) {
                unsigned char aip = (unsigned char)a[i-1];
                if ((int)ai + 1 > (int)aip) continue;
            }
            // 新单元 (i, j=ai)
            int j = (int)ai;
            int w = weight(i, j, m, h);
            long long step_cost = cfg_.dV * (long long)w * (long long)(level + 1);

            PackedFrontier nxt = key;
            nxt.bytes[i] = (char)(ai + 1);

            int pid = idmap.count(key) ? idmap[key] : -1;
            push_state(nxt, cdist + step_cost, level + 1, pid, i, j);
        }

        ++expanded;
        if (cfg_.verbose && (expanded % 500000 == 0)) {
            cout << "  [progress] expanded=" << expanded
                 << " pq=" << pq.size()
                 << " cached=" << best.size() << endl;
        }
    }

    if (Z < 0) {
        throw std::runtime_error("FrontierMcmfColumnPlacer: no path found (unexpected).");
    }

    // ---------- 回溯构造排名矩阵 ----------
    vector<vector<int>> Y(m, vector<int>(h, 0));
    vector<pair<int,int>> seq; seq.reserve(N);

    PackedFrontier curkey = goal;
    while (true) {
        const NodeRec& rec = best[curkey];
        if (rec.put_i >= 0) {
            seq.emplace_back(rec.put_i, rec.put_j);
        }
        if (rec.parent_id < 0) break;
        curkey = id2key[rec.parent_id];
    }
    // seq 目前从最后一步到第一步，倒转并填 rank
    std::reverse(seq.begin(), seq.end());
    for (int t = 0; t < (int)seq.size(); ++t) {
        int i = seq[t].first;
        int j = seq[t].second;
        Y[i][j] = t + 1;
    }

    FrontierResult R;
    R.m = m; R.h = h; R.dV = cfg_.dV;
    R.y_order = std::move(Y);
    R.x_of_col.assign(h, 0);
    R.mcmf_cost = Z;
    R.actual_hpwl = hpwl_edges_sum(R.y_order, R.x_of_col, cfg_.dV);
    R.oc_ok = check_global_OC(R.y_order);

    // 唯一占用 (x,y)
    {
        std::unordered_set<uint64_t> occ; occ.reserve((size_t)N*2);
        bool ok = true;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < h; ++j) {
                uint64_t x = 0ull;
                uint64_t yv = (uint64_t)R.y_order[i][j];
                uint64_t key64 = (x << 32) ^ yv;
                if (!occ.insert(key64).second) ok = false;
            }
        R.unique_ok = ok;
    }

    R.states_expanded = expanded;
    R.states_cached   = best.size();

    if (cfg_.verbose) {
        cout << "[Pure-MCMF] states_expanded=" << R.states_expanded
             << ", states_cached=" << R.states_cached << "\n";
        cout << "  [verify] MCMF cost=" << R.mcmf_cost
             << ", Actual HPWL=" << R.actual_hpwl
             << ", OC=" << (R.oc_ok ? "OK" : "FAIL")
             << ", Unique=" << (R.unique_ok ? "OK" : "FAIL") << endl;
    }

    return R;
}
