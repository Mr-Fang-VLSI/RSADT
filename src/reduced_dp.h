#ifndef REDUCED_DP_H
#define REDUCED_DP_H

#include <vector>
#include <string>
#include <unordered_map>
#include <limits>
#include <cstdint>

struct ReducedDPConfig {
    bool verbose = false;          // 中间打印
    bool print_layout = true;      // 打印最终布局
    bool enable_left_chain = true; // 保守：左列开行串（-1 包）
    bool enable_right_chain = true;// 保守：右列收尾串（+1 包）
    bool canonize_start = true;    // 起点也做 canonicalization
    bool macro_col_in_canon = true;// canonicalization 里是否做列级级联
    int  max_col_cascade  = 1;     // 每次 canonicalize 最多吃几列（安全：默认 1）
};

struct ReducedDPResult {
    long long total_cost = 0;
    std::vector<std::vector<int>> y_order; // m x h, 1..N
    size_t num_states = 0; // 归并后的节点数
    size_t num_edges = 0;
    double states_over_m2h2 = 0.0;
    double states_over_binom = 0.0; // 相对 C(m+h, m)
};

class ReducedDP {
public:
    ReducedDP(int m, int h, const ReducedDPConfig& cfg);

    ReducedDPResult solve();

private:
    struct State {
        std::vector<uint16_t> a; // 每行已放置数（0..h），非增
        bool operator==(const State& o) const { return a == o.a; }
    };
    struct StateHasher {
        size_t operator()(const State& s) const {
            uint64_t h = 1469598103934665603ULL;
            for(uint16_t x : s.a){
                h ^= (uint64_t)x + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2);
            }
            return (size_t)h;
        }
    };
    struct EdgeInfo {
        enum Type : uint8_t { STEP=0, MACRO_COL=1, MACRO_LEFT=2, MACRO_RIGHT=3, COMPOSITE=4 } type;
        int row=-1, col=-1;  // STEP
        int r= -1;           // MACRO_COL 的列号
        int k= 0;            // MACRO_LEFT/RIGHT 的条数
        int seq_id = -1;     // COMPOSITE：序列存储索引
    };

    // 基本尺寸
    int m_, h_, N_;
    ReducedDPConfig cfg_;
    // phi(i,j) : 等权
    std::vector<std::vector<int>> phi_;

    // 存放复合边的展开序列（用于回溯）
    std::vector<std::vector<EdgeInfo>> seq_store_;

    // ===== 基本工具 =====
    void build_phi();
    static long long hpwl_equal(const std::vector<std::vector<int>>& y);
    static bool check_oc(const std::vector<std::vector<int>>& y);
    static inline int level_of(const State& s){ int L=0; for(uint16_t x: s.a) L+= (int)x; return L; }
    static inline bool step_feasible(const State& s, int i){
        if(s.a[i] >= std::numeric_limits<uint16_t>::max()) return false;
        if(i==0) return true;
        return s.a[i] < s.a[i-1];
    }

    // ===== 无损宏规则（单列/保守串） =====
    bool make_macro_col(const State& s, int L0, long long& cost, State& t, EdgeInfo& einfo) const;
    bool make_macro_left(const State& s, int L0, long long& cost, State& t, EdgeInfo& einfo) const;
    bool make_macro_right(const State& s, int L0, long long& cost, State& t, EdgeInfo& einfo) const;

    // 单步 exact 兜底
    void expand_unit_steps(const State& s, int L0,
                           std::vector<std::pair<State,long long>>& outs,
                           std::vector<EdgeInfo>& eouts) const;

    // ===== Canonicalization：最大化“安全宏推进”，把整段状态等价归并到代表状态 =====
    void canonicalize(State& s, int L0, long long& delta, std::vector<EdgeInfo>& seq) const;

    // 初始/终止状态
    State start_state() const;
    State goal_state() const;

    // 回溯：把（可能是 COMPOSITE 的）边展开为坐标序列，生成 y
    void reconstruct_y(const std::vector<EdgeInfo>& path_edges, std::vector<std::vector<int>>& y) const;

    static void print_layout(const std::vector<std::vector<int>>& y);

    // 组合数 C(n,k) 的近似（double），用于 states_over_binom
    static double binom_approx(int n, int k);
};

#endif // REDUCED_DP_H
