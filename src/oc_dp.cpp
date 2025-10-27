#include "oc_dp.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>

// ===== 工具 =====
int OCDP::ceil_log2(int x){
    int k=0; int v=1;
    while(v < x){ v<<=1; ++k; }
    return k;
}

// ===== HPWL / OC 校验 =====
long long OCDP::hpwl_sum_equal(const std::vector<std::vector<int>>& y){
    int m=(int)y.size(), h=(int)y[0].size();
    long long W=0;
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        if(i+1<m) W += llabs((long long)y[i+1][j]-y[i][j]);
        if(j+1<h) W += llabs((long long)y[i][j+1]-y[i][j]);
    }
    return W;
}
long long OCDP::hpwl_sum_weighted(const std::vector<std::vector<int>>& y, const OCWeights& W){
    int m=(int)y.size(), h=(int)y[0].size();
    long long S=0;
    if(!W.enabled) return hpwl_sum_equal(y);
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        if(i+1<m) S += llabs((long long)y[i+1][j]-y[i][j]) * W.wV[i][j];
        if(j+1<h) S += llabs((long long)y[i][j+1]-y[i][j]) * W.wH[i][j];
    }
    return S;
}
bool OCDP::check_OC(const std::vector<std::vector<int>>& y){
    int m=(int)y.size(), h=(int)y[0].size(), N=m*h;
    std::vector<std::pair<int,int>> pos(N+1, {-1,-1});
    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        int r = y[i][j];
        if(r<1 || r>N) return false;
        pos[r] = {i,j};
    }
    std::vector<std::vector<char>> seen(m, std::vector<char>(h, 0));
    for(int r=1;r<=N;++r){
        auto [i,j]=pos[r];
        if(i>0 && !seen[i-1][j]) return false;
        if(j>0 && !seen[i][j-1]) return false;
        seen[i][j]=1;
    }
    return true;
}

// ===== phi 构建：phi = (上+左) - (下+右) =====
void OCDP::build_phi(const OCWeights& W){
    phi_.assign(m_, std::vector<long long>(h_, 0));
    if(!W.enabled){
        for(int i=0;i<m_;++i){
            for(int j=0;j<h_;++j){
                long long up    = (i>0)?1:0;
                long long left  = (j>0)?1:0;
                long long down  = (i<m_-1)?1:0;
                long long right = (j<h_-1)?1:0;
                phi_[i][j] = (up+left) - (down+right);
            }
        }
    }else{
        for(int i=0;i<m_;++i){
            for(int j=0;j<h_;++j){
                long long up    = (i>0    ? W.wV[i-1][j] : 0);
                long long left  = (j>0    ? W.wH[i][j-1] : 0);
                long long down  = (i<m_-1 ? W.wV[i][j]   : 0);
                long long right = (j<h_-1 ? W.wH[i][j]   : 0);
                phi_[i][j] = (up+left) - (down+right);
            }
        }
    }
}

// ===== 打包初始化 =====
void OCDP::init_packers(){
    SEG_  = cfg_.use_pow2_pack ? ceil_log2(h_+1) : 0;
    if(cfg_.use_pow2_pack){
        MASK_ = (SEG_==64? ~0ULL : ((1ULL<<SEG_)-1ULL));
        BITS_ = m_ * SEG_;
    }else{
        BASE_ = (long long)(h_+1);
        BITS_ = 0; // 用不到
    }
    nlimbs_ = (cfg_.use_pow2_pack ? ((BITS_+63)/64) : 4); // radix 用 4 limbs
    if(nlimbs_<=0) nlimbs_=1;
    if(nlimbs_>4)  nlimbs_=4;

    off_bit_.assign(m_, 0);
    if(cfg_.use_pow2_pack){
        for(int i=0;i<m_;++i) off_bit_[i] = i*SEG_;
    }else{
        // 预计算 powB
        for(int l=0;l<4;++l) powB_[l].assign(m_, 0);
        // powB[idx] = BASE^idx  （256-bit，小端）
        auto mul_small = [&](uint64_t a[4], long long B){
            __uint128_t carry = 0;
            for(int p=0;p<4;++p){
                __uint128_t v = (__uint128_t)a[p]* (unsigned long long)B + carry;
                a[p] = (uint64_t)(v & 0xFFFFFFFFFFFFFFFFULL);
                carry = (v >> 64);
            }
        };
        uint64_t cur[4]={1,0,0,0};
        for(int idx=0; idx<m_; ++idx){
            for(int p=0;p<4;++p) powB_[p][idx] = cur[p];
            mul_small(cur, BASE_);
        }
    }
}

// ===== pow2 读写 =====
int OCDP::get_digit_pow2(const KeyX& k, int idx) const{
    int ob = off_bit_[idx];
    int li = ob >> 6;
    int sh = ob & 63;
    if(sh + SEG_ <= 64){
        uint64_t w = k.v[li];
        return (int)((w >> sh) & MASK_);
    }else{
        int r = (sh + SEG_) - 64;
        uint64_t lo = (k.v[li] >> sh);
        uint64_t hi = (li+1<4 ? k.v[li+1] : 0ULL) & ((1ULL<<r)-1ULL);
        return (int)((lo | (hi << (64-sh))) & MASK_);
    }
}
void OCDP::set_digit_pow2(KeyX& k, int idx, int val) const{
    int ob = off_bit_[idx];
    int li = ob >> 6;
    int sh = ob & 63;
    uint64_t m = (SEG_==64? ~0ULL : ((1ULL<<SEG_)-1ULL));
    if(sh + SEG_ <= 64){
        uint64_t clr = ~(m << sh);
        k.v[li] = (k.v[li] & clr) | ((uint64_t)val << sh);
    }else{
        int r = (sh + SEG_) - 64;
        // low part
        uint64_t low_mask = ~((~0ULL) << sh);
        k.v[li] = (k.v[li] & low_mask) | ((uint64_t)val << sh);
        // high part
        uint64_t high_mask = ~((1ULL<<r)-1ULL);
        uint64_t hi = (uint64_t)val >> (64 - sh);
        if(li+1<4) k.v[li+1] = (k.v[li+1] & high_mask) | (hi & ((1ULL<<r)-1ULL));
    }
}
OCDP::KeyX OCDP::inc_digit_pow2(const KeyX& k, int idx) const{
    KeyX t = k;
    int v = get_digit_pow2(k, idx);
    set_digit_pow2(t, idx, v+1);
    return t;
}
OCDP::KeyX OCDP::dec_digit_pow2(const KeyX& k, int idx) const{
    KeyX t = k;
    int v = get_digit_pow2(k, idx);
    set_digit_pow2(t, idx, v-1);
    return t;
}

// ===== radix 读写 =====
void OCDP::add_big(KeyX& out, const KeyX& a, const std::vector<uint64_t> add[4], int idx) const{
    __uint128_t carry=0;
    for(int p=0;p<nlimbs_;++p){
        __uint128_t v = (__uint128_t)a.v[p] + add[p][idx] + carry;
        out.v[p] = (uint64_t)v;
        carry = (v>>64);
    }
    for(int p=nlimbs_; p<4; ++p) out.v[p]=0;
}
OCDP::KeyX OCDP::inc_digit_radix(const KeyX& k, int idx) const{
    KeyX t; add_big(t, k, powB_, idx); return t;
}
OCDP::KeyX OCDP::dec_digit_radix(const KeyX& k, int idx) const{
    KeyX t{{0,0,0,0}};
    __int128 carry=0;
    for(int p=0;p<nlimbs_;++p){
        __int128 v = (__int128)k.v[p] - ( __int128)powB_[p][idx] + carry;
        t.v[p] = (uint64_t)v;
        carry = (v<0? -1:0);
    }
    for(int p=nlimbs_; p<4; ++p) t.v[p]=0;
    return t;
}
int OCDP::get_digit_radix_slow(const KeyX& k, int idx) const{
    // 计算 ((k / BASE^idx) % BASE)
    // 复制 k
    uint64_t num[4] = {k.v[0],k.v[1],k.v[2],k.v[3]};
    // 连续 idx 次：num = num / BASE
    for(int t=0;t<idx;++t){
        __uint128_t rem=0;
        for(int p=3;p>=0;--p){
            __uint128_t cur = (rem<<64) | num[p];
            uint64_t q = (uint64_t)(cur / (unsigned long long)BASE_);
            uint64_t r = (uint64_t)(cur % (unsigned long long)BASE_);
            num[p] = q; rem = r;
        }
    }
    // 现在取 num % BASE
    __uint128_t rem2=0;
    for(int p=3;p>=0;--p){
        __uint128_t cur = (rem2<<64) | num[p];
        rem2 = cur % (unsigned long long)BASE_;
    }
    return (int)rem2;
}

// ===== Key 比较 =====
int OCDP::cmp_key(const KeyX& a, const KeyX& b) const{
    for(int p=nlimbs_-1; p>=0; --p){
        if(a.v[p] < b.v[p]) return -1;
        if(a.v[p] > b.v[p]) return  1;
    }
    return 0;
}

// ===== 前向推进到 L0+steps =====
void OCDP::forward_to_level(const KeyX& key0, int L0, int steps, std::vector<KV>& out) const{
    std::vector<KV> cur, nxt;
    cur.reserve(1024);
    cur.push_back({key0, 0});

    for(int lvl=0; lvl<steps; ++lvl){
        nxt.clear();
        nxt.reserve(cur.size() * (size_t)std::min(m_, h_) + 16);
        for(const auto& kv : cur){
            int ai[64];
            for(int i=0;i<m_;++i) ai[i]=get_digit(kv.key,i);
            for(int i=0;i<m_;++i){
                if(ai[i] >= h_) continue;
                if(i>0 && ai[i]+1 > ai[i-1]) continue;
                long long c = (long long)(L0 + lvl + 1) * phi_[i][ai[i]];
                KeyX k2 = inc_digit(kv.key, i);
                nxt.push_back({k2, kv.dist + c});
            }
        }
        std::sort(nxt.begin(), nxt.end(), CmpKV{this});
        std::vector<KV> ded; ded.reserve(nxt.size());
        for(size_t i=0;i<nxt.size();){
            size_t j=i+1;
            long long bestd = nxt[i].dist;
            KeyX bestk = nxt[i].key;
            while(j<nxt.size() && cmp_key(nxt[j].key, nxt[i].key)==0){
                if(nxt[j].dist < bestd){ bestd = nxt[j].dist; bestk = nxt[j].key; }
                ++j;
            }
            ded.push_back({bestk, bestd});
            i=j;
        }
        cur.swap(ded);
        if(cfg_.progress && ((L0+lvl+1)%16==0))
            std::cerr<<"  [F] reach L="<<(L0+lvl+1)<<" states="<<cur.size()<<"\n";
    }
    out.swap(cur);
}

// ===== 后向推进到 Lhi-steps =====
void OCDP::backward_to_level(const KeyX& keyN, int LN, int steps, std::vector<KV>& out) const{
    std::vector<KV> cur, nxt;
    cur.reserve(1024);
    cur.push_back({keyN, 0});

    for(int t=0; t<steps; ++t){
        int Lcur = LN - t;
        nxt.clear();
        nxt.reserve(cur.size() * (size_t)std::min(m_, h_) + 16);
        for(const auto& kv : cur){
            int ai[64];
            for(int i=0;i<m_;++i) ai[i]=get_digit(kv.key,i);
            for(int i=0;i<m_;++i){
                if(ai[i] <= 0) continue;
                if(i<m_-1 && (ai[i]-1) < ai[i+1]) continue;
                int jprev = ai[i]-1;
                long long c = (long long)Lcur * phi_[i][jprev];
                KeyX k2 = dec_digit(kv.key, i);
                nxt.push_back({k2, kv.dist + c});
            }
        }
        std::sort(nxt.begin(), nxt.end(), CmpKV{this});
        std::vector<KV> ded; ded.reserve(nxt.size());
        for(size_t i=0;i<nxt.size();){
            size_t j=i+1;
            long long bestd = nxt[i].dist;
            KeyX bestk = nxt[i].key;
            while(j<nxt.size() && cmp_key(nxt[j].key, nxt[i].key)==0){
                if(nxt[j].dist < bestd){ bestd = nxt[j].dist; bestk = nxt[j].key; }
                ++j;
            }
            ded.push_back({bestk, bestd});
            i=j;
        }
        cur.swap(ded);
        if(cfg_.progress && ((Lcur-1)%16==0))
            std::cerr<<"  [B] reach L="<<(Lcur-1)<<" states="<<cur.size()<<"\n";
    }
    out.swap(cur);
}

// ===== 小区间直接求 & 写入 y（强制终点 = keyT） =====
void OCDP::solve_small_and_fill(const KeyX& key0, const KeyX& keyT,
                                int L0, int steps,
                                std::vector<std::vector<int>>& y) const{
    struct Node { KeyX key; long long dist; };
    std::vector<std::vector<Node>> layers(steps+1);
    std::vector<std::vector<int>>  parent_idx(steps+1);
    std::vector<std::vector<int>>  parent_row(steps+1);

    layers[0].push_back({key0, 0});

    for(int s=0; s<steps; ++s){
        // 生成候选（key,dist,parentIdx,row）
        std::vector<std::tuple<KeyX,long long,int,int>> cand;
        cand.reserve(layers[s].size() * (size_t)std::min(m_, h_) + 8);

        for(size_t pid=0; pid<layers[s].size(); ++pid){
            const auto& n = layers[s][pid];
            int ai[64]; for(int i=0;i<m_;++i) ai[i]=get_digit(n.key,i);
            for(int i=0;i<m_;++i){
                if(ai[i]>=h_) continue;
                if(i>0 && ai[i]+1>ai[i-1]) continue;
                long long c = (long long)(L0 + s + 1) * phi_[i][ai[i]];
                KeyX k2 = inc_digit(n.key, i);
                cand.emplace_back(k2, n.dist + c, (int)pid, i);
            }
        }
        if(cand.empty()){
            if(cfg_.verbose)
                std::cerr<<"[smallDP] empty candidate at s="<<s<<", steps="<<steps<<"\n";
            return;
        }
        // sort + unique(min) 并且保留父信息
        std::sort(cand.begin(), cand.end(), [&](auto& A, auto& B){
            const KeyX& x = std::get<0>(A);
            const KeyX& yk= std::get<0>(B);
            int c = cmp_key(x, yk);
            if(c!=0) return c<0;
            return std::get<1>(A) < std::get<1>(B);
        });
        layers[s+1].clear(); parent_idx[s+1].clear(); parent_row[s+1].clear();
        layers[s+1].reserve(cand.size());
        for(size_t i=0;i<cand.size();){
            size_t j=i+1;
            auto best = cand[i];
            while(j<cand.size() && cmp_key(std::get<0>(cand[j]), std::get<0>(cand[i]))==0){
                if(std::get<1>(cand[j]) < std::get<1>(best)) best=cand[j];
                ++j;
            }
            layers[s+1].push_back({ std::get<0>(best), std::get<1>(best) });
            parent_idx[s+1].push_back( std::get<2>(best) );
            parent_row[s+1].push_back( std::get<3>(best) );
            i=j;
        }
        if(cfg_.verbose && ((L0+s+1)%16==0))
            std::cerr<<"  [smallDP] L="<<(L0+s+1)<<" states="<<layers[s+1].size()<<"\n";
    }

    // 末层：优先找 key == keyT
    int cur = -1;
    for(size_t i=0;i<layers[steps].size();++i){
        if(cmp_key(layers[steps][i].key, keyT)==0){ cur=(int)i; break; }
    }
    if(cur<0){
        // 理论上不该发生；回退为选最小 dist 并告警
        if(cfg_.verbose){
            std::cerr<<"[smallDP][warn] target key not found in last layer; fallback to min-dist node.\n";
        }
        long long bestd = layers[steps][0].dist; cur=0;
        for(size_t i=1;i<layers[steps].size();++i)
            if(layers[steps][i].dist < bestd){ bestd = layers[steps][i].dist; cur=(int)i; }
    }

    // 回溯 rows
    std::vector<int> rows; rows.reserve(steps);
    for(int s=steps; s>=1; --s){
        int r = parent_row[s][cur];
        rows.push_back(r);
        cur = parent_idx[s][cur];
    }
    std::reverse(rows.begin(), rows.end());

    // 重放以确定 (i,j) 并写入 y
    int ai[64]; for(int i=0;i<m_;++i) ai[i]=get_digit(key0,i);
    int rank = L0;
    for(int t=0;t<steps;++t){
        int i = rows[t];
        int j = ai[i];
        y[i][j] = ++rank;
        ai[i]++;
    }
}

// ===== Hirschberg =====
void OCDP::reconstruct(int Llo, int Lhi, const KeyX& key_lo, const KeyX& key_hi,
                       std::vector<std::vector<int>>& y) const{
    int len = Lhi - Llo;
    if(len<=0) return;
    if(len <= cfg_.small_thresh){
        // **修正点**：小区间必须从 key_lo 到 key_hi
        solve_small_and_fill(key_lo, key_hi, Llo, len, y);
        return;
    }
    int Lmid = Llo + len/2;
    int stepsF = Lmid - Llo;
    int stepsB = Lhi - Lmid;

    std::vector<KV> F, B;
    forward_to_level(key_lo, Llo, stepsF, F);
    backward_to_level(key_hi, Lhi, stepsB, B);

    long long best = (1LL<<62);
    KeyX bestK{{0,0,0,0}};
    size_t i=0, j=0;
    while(i<F.size() && j<B.size()){
        int c = cmp_key(F[i].key, B[j].key);
        if(c==0){
            long long v = F[i].dist + B[j].dist;
            if(v < best){ best=v; bestK = F[i].key; }
            ++i; ++j;
        }else if(c<0) ++i;
        else ++j;
    }
    if(best==(1LL<<62)){
        // 不应发生；输出诊断后直接小区间求解兜底
        if(cfg_.verbose){
            std::cerr<<"[reconstruct][warn] no matching mid key; fallback small solve len="<<len<<"\n";
        }
        solve_small_and_fill(key_lo, key_hi, Llo, len, y);
        return;
    }

    // 左半 & 右半递归
    reconstruct(Llo, Lmid, key_lo, bestK, y);
    reconstruct(Lmid, Lhi, bestK, key_hi, y);
}

// ===== solve =====
OCResult OCDP::solve(int m, int h, const OCWeights& W){
    m_=m; h_=h; N_=m*h;
    build_phi(W);
    init_packers();

    // 起止 key
    KeyX key0 = key_zero();
    KeyX keyN = key_zero();
    for(int i=0;i<m_;++i){
        if(cfg_.use_pow2_pack) set_digit_pow2(keyN, i, h_);
        else{
            for(int t=0;t<h_;++t) keyN = inc_digit_radix(keyN, i);
        }
    }

    std::vector<std::vector<int>> y(m_, std::vector<int>(h_,0));

    reconstruct(0, N_, key0, keyN, y);

    // 计算总代价
    long long total=0;
    int L=0;
    std::vector<std::pair<int,int>> pos(N_+1,{0,0});
    for(int i=0;i<m_;++i) for(int j=0;j<h_;++j) pos[y[i][j]]={i,j};
    for(int r=1;r<=N_;++r){
        auto [i,j]=pos[r];
        total += (long long)(L+1) * phi_[i][j];
        ++L;
    }

    OCResult R;
    R.total_cost = total;
    R.y_order = std::move(y);
    return R;
}
