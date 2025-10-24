#pragma once
#include <vector>
#include <cstdint>
#include <array>
#include <unordered_map>
#include <stdexcept>
#include <limits>

struct OCShortestResult {
    int m = 0, h = 0, n = 0;
    long long dV = 1;
    std::vector<std::vector<int>> y_order;
    long long total_cost = 0; // = HPWL (dV=1) for verification
};

struct Weights {
    // H: m x (h-1), V: (m-1) x h
    std::vector<std::vector<long long>> wH, wV;
    bool enabled = false;
};

class lightOCShortest {
public:
    struct Config {
        long long dV = 1;
        bool verbose = true;
        int  max_steps = -1;     // -1 -> full (m*h)
        bool progress  = false;
        int  omp_threads = 0;    // 0 -> all cores
        bool low_mem = true;     // Hirschberg low-memory rec.
        bool use_pow2_pack = true;
        Weights W;
    };
    explicit lightOCShortest(const Config& cfg) : cfg_(cfg) {}
    OCShortestResult solve(int m, int h);

    // ==== public helpers (verification) ====
    static long long hpwl_sum_equal(const std::vector<std::vector<int>>& y, long long dV);
    static long long hpwl_sum_weighted(const std::vector<std::vector<int>>& y,
                                       const Weights& W, long long dV);
    static bool check_OC_lin(const std::vector<std::vector<int>>& y);

    // ==== 64-bit base-(h+1) digit ops (for small cases) ====
    static inline uint64_t encode_digit_inc64(uint64_t key, int idx, const std::vector<uint64_t>& powB) {
        return key + powB[idx];
    }
    static inline uint64_t encode_digit_dec64(uint64_t key, int idx, const std::vector<uint64_t>& powB) {
        return key - powB[idx];
    }
    static inline int digit_at64(uint64_t key, int idx, const std::vector<uint64_t>& powB, uint64_t B) {
        return int((key / powB[idx]) % B);
    }

    // ==== 2^k bit-packed key (supports up to 256 bits) ====
    struct KeyX {
        std::array<uint64_t,4> w{}; // only first nwords_ used
        bool operator==(const KeyX& o) const {
            return w[0]==o.w[0] && w[1]==o.w[1] && w[2]==o.w[2] && w[3]==o.w[3];
        }
    };
    struct KeyXHasher {
        size_t operator()(const KeyX& k) const noexcept {
            uint64_t x = k.w[0] ^ (k.w[1]*0x9e3779b185ebca87ULL)
                       ^ (k.w[2]<<1) ^ (k.w[3]*0xbf58476d1ce4e5b9ULL);
            return (size_t)(x ^ (x>>33));
        }
    };

    // bit-pack accessors (defined in .cpp)
    int  get_digit(const KeyX& k, int idx) const;
    void set_digit(KeyX& k, int idx, int v) const;
    KeyX inc_digit(const KeyX& k, int idx) const;
    KeyX dec_digit(const KeyX& k, int idx) const;

    // expose for free helpers
    inline long long phi_at(int i,int j) const { return phi_[i][j]; }

    // reconstruct (declared here, implemented in .cpp)
    void reconstruct_hirschberg64(int m,int h,long long dV,
        const std::vector<uint64_t>& powB,uint64_t B,bool progress,int omp_threads,
        int l_lo,uint64_t key_lo,int l_hi,uint64_t key_hi,std::vector<std::vector<int>>& y);

    void reconstruct_hirschbergX(int m,int h,long long dV,
        bool progress,int omp_threads,
        int l_lo,KeyX key_lo,int l_hi,KeyX key_hi,std::vector<std::vector<int>>& y);

private:
    Config cfg_;

    // phi precompute (m x h)
    std::vector<std::vector<long long>> phi_;
    void build_phi(int m,int h);

    // bit-pack codec params
    int  seg_ = 0;                 // SEG = ceil(log2(h+1))
    int  total_bits_ = 0;          // m * SEG
    int  nwords_ = 1;              // 1..4
    uint64_t mask_ = 0;            // (1<<SEG)-1 (SEG=64 handled in .cpp)
    std::array<int, 32> shift_{};  // SHIFT[i] = (m-1-i)*SEG
    void init_bitpack_codec(int m,int h);
};
