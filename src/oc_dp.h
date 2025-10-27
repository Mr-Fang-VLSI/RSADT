#ifndef OC_DP_H
#define OC_DP_H

#include <vector>
#include <cstdint>
#include <string>
#include <utility>

struct OCSolveConfig {
    bool use_pow2_pack = true;   // 2^k 位宽打包（默认）。若关，则用混合基数 (h+1) 打包
    bool progress       = true;  // 分层进度
    bool verbose        = false; // 详细日志
    int  small_thresh   = 64;    // 小区间长度阈值（≤此值走“局部父指针 DP”，强制终点为 key_hi）
};

struct OCWeights {
    bool enabled = false; // false => 等权
    // wH[i][j] : (i,j) -> (i,j+1)  大小 m x (h-1)
    // wV[i][j] : (i,j) -> (i+1,j)  大小 (m-1) x h
    std::vector<std::vector<long long>> wH, wV;
};

struct OCResult {
    long long total_cost = 0;
    std::vector<std::vector<int>> y_order; // m x h, 1..N
};

class OCDP {
public:
    explicit OCDP(const OCSolveConfig& cfg): cfg_(cfg) {}

    // 主入口
    OCResult solve(int m, int h, const OCWeights& W);

    // 工具：验证
    static long long hpwl_sum_equal(const std::vector<std::vector<int>>& y);
    static long long hpwl_sum_weighted(const std::vector<std::vector<int>>& y, const OCWeights& W);
    static bool check_OC(const std::vector<std::vector<int>>& y);

private:
    // ===== 关键：可变位宽 Key 打包（最多 256bit, 4 limbs）=====
    struct KeyX {
        uint64_t v[4]; // 按 nlimbs_ 使用
        bool operator==(const KeyX& o) const {
            return v[0]==o.v[0] && v[1]==o.v[1] && v[2]==o.v[2] && v[3]==o.v[3];
        }
    };
    struct KV { KeyX key; long long dist; };

    OCSolveConfig cfg_;
    int m_=0, h_=0, N_=0;

    // phi(i,j) = (上+左) - (下+右)，边界越界视 0；加权直接把权代入
    std::vector<std::vector<long long>> phi_;

    // 打包参数
    int SEG_=0;                 // 每槽位宽（bit），pow2 模式：ceil_log2(h+1)
    int BITS_=0;                // 总位宽 m*SEG
    int nlimbs_=0;              // 1..4
    uint64_t MASK_=0;           // (1<<SEG)-1
    // 对 pow2 模式的位偏移
    std::vector<int> off_bit_;  // 每槽 bit offset
    // 对混合基数 (h+1) 模式
    long long BASE_=0;
    std::vector<uint64_t> powB_[4]; // powB_[limb][idx]，按 nlimbs_ 使用

    // ===== KeyX 工具 =====
    KeyX key_zero() const { KeyX k{{0,0,0,0}}; return k; }
    void init_packers();

    // pow2 读写
    int  get_digit_pow2 (const KeyX& k, int idx) const;
    void set_digit_pow2 (KeyX& k, int idx, int v) const;
    KeyX inc_digit_pow2 (const KeyX& k, int idx) const;
    KeyX dec_digit_pow2 (const KeyX& k, int idx) const;

    // radix 读写（基数 = h+1）
    void add_big(KeyX& out, const KeyX& a, const std::vector<uint64_t> add[4], int idx) const;
    KeyX inc_digit_radix(const KeyX& k, int idx) const;
    KeyX dec_digit_radix(const KeyX& k, int idx) const;
    int  get_digit_radix_slow(const KeyX& k, int idx) const; // 仅备用
    int  get_digit(const KeyX& k, int idx) const {
        return cfg_.use_pow2_pack ? get_digit_pow2(k,idx) : get_digit_radix_slow(k,idx);
    }
    KeyX inc_digit(const KeyX& k, int idx) const {
        return cfg_.use_pow2_pack ? inc_digit_pow2(k,idx) : inc_digit_radix(k,idx);
    }
    KeyX dec_digit(const KeyX& k, int idx) const {
        return cfg_.use_pow2_pack ? dec_digit_pow2(k,idx) : dec_digit_radix(k,idx);
    }

    // 比较器/相等（按 nlimbs_）
    int cmp_key(const KeyX& a, const KeyX& b) const;
    struct CmpKV {
        const OCDP* S;
        bool operator()(const KV& x, const KV& y) const { return S->cmp_key(x.key,y.key)<0; }
    };

    // ===== DP 层推进（前向/后向）=====
    void forward_to_level(const KeyX& key0, int L0, int steps, std::vector<KV>& out) const;
    void backward_to_level(const KeyX& keyN, int LN, int steps, std::vector<KV>& out) const;

    // 小区间直接重建（父指针）并写入 y —— **强制终点等于 keyT**
    void solve_small_and_fill(const KeyX& key0, const KeyX& keyT,
                              int L0, int steps,
                              std::vector<std::vector<int>>& y) const;

    // Hirschberg 二分重建
    void reconstruct(int Llo, int Lhi, const KeyX& key_lo, const KeyX& key_hi,
                     std::vector<std::vector<int>>& y) const;

    // 构建 phi
    void build_phi(const OCWeights& W);

    // 工具
    static int ceil_log2(int x);
};

#endif // OC_DP_H
