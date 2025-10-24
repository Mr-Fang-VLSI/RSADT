#pragma once
#include <vector>
#include <string>
#include <cstdint>

struct OCClosureResult {
    int m=0, h=0, n=0;
    long long dV=1;
    std::vector<std::vector<int>> y_order; // rank matrix [m][h], values in [1..n]
    long long total_cost=0;                 // = HPWL (验证等价于 sum phi*rank)
};

struct OCWeights {
    // H: m x (h-1), V: (m-1) x h
    std::vector<std::vector<long long>> wH, wV;
    bool enabled=false;
};

class OCClosure {
public:
    struct Config {
        long long dV = 1;          // 代价缩放
        bool verbose = true;        // 关键步骤日志
        bool progress = false;      // Dinkelbach 迭代日志
    };

    explicit OCClosure(const Config& cfg): cfg_(cfg) {}

    // 输入：m,h（网格），可选权重 W（不传或 enabled=false 表示等权）
    // 输出：全局最优 y（严格满足 OC），total_cost 为 HPWL（加权/等权）
    OCClosureResult solve(int m, int h, const OCWeights& W);

    // 验证与工具
    static bool check_OC_lin(const std::vector<std::vector<int>>& y);
    static long long hpwl_sum_equal(const std::vector<std::vector<int>>& y, long long dV);
    static long long hpwl_sum_weighted(const std::vector<std::vector<int>>& y,
                                       const OCWeights& W, long long dV);

private:
    Config cfg_;

    // ϕ(i,j) = (up + left) - (down + right) ；越界边=0；乘以 dV
    std::vector<std::vector<long long>> phi_;
    void build_phi(int m,int h, const OCWeights& W);

    // ---- Lawler–Sidney 分解(每轮找“最大密度”的向上闭合集，并排在最后) ----
    // 整体递归：对 mask（子图）返回“从早到晚”的节点顺序 id 列表
    std::vector<int> solve_recursive(const std::vector<char>& mask) const;

    // Dinkelbach：最大化  A(S)/|S| ，S 为“向上闭合（upward‑closed）”且非空
    // 返回选中的 S（用 maskS 标出），并返回 A=∑w，B=|S|
    void find_max_density_upward_closed(const std::vector<char>& mask,
                                        std::vector<char>& maskS,
                                        long long& A_sum, int& B_count) const;

    // ---- 最大权闭包：max ∑ q_i , S upward-closed ----
    // q_i 由 Dinkelbach 的 (w_i * B - A) 或字典序放大后的 q'_i 得到（64位，内部用 __int128 防溢出）
    // 返回 closure S（maskS）与 sumQ
    long long max_weight_upward_closure(const std::vector<char>& mask,
                                        const std::vector<long long>& q,
                                        std::vector<char>& maskS) const;

    // ---- 图结构（原始 DAG：每个点连向“下/右”） ----
    int m_=0, h_=0, N_=0;
    inline int id(int i,int j) const { return i*h_ + j; }
    std::vector<std::vector<int>> succ_; // upward-closure 约束边：u -> v（右/下）

    // ---- 内部工具 ----
    static void complement_mask(const std::vector<char>& mask,
                                const std::vector<char>& S,
                                std::vector<char>& maskLeft);

    long long sum_w_on_mask(const std::vector<char>& mask) const;
    int count_on_mask(const std::vector<char>& mask) const;

    // 取 mask 上的 “sink（无后继）集合”，必为 upward-closed，且非空（只要 mask 非空）
    void sinks_upward_closed(const std::vector<char>& mask,
                             std::vector<char>& S_sink) const;

    // 调试：统计 S 是否 upward-closed（仅 verbose 时用）
    bool is_upward_closed_on_mask(const std::vector<char>& mask,
                                  const std::vector<char>& S) const;
};
