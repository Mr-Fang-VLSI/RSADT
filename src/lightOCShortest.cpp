#include "lightOCShortest.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <tuple>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

using std::vector;
using std::unordered_map;
using std::cout;
using std::endl;

// ===== helpers: verification =====
long long lightOCShortest::hpwl_sum_equal(const vector<vector<int>>& y, long long dV){
    const int m=(int)y.size(), h=(int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j)
        S += (long long)std::llabs((long long)y[i][j+1]-y[i][j]) * dV;
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j)
        S += (long long)std::llabs((long long)y[i+1][j]-y[i][j]) * dV;
    return S;
}
long long lightOCShortest::hpwl_sum_weighted(const vector<vector<int>>& y,
                                             const Weights& W, long long dV){
    if(!W.enabled) return hpwl_sum_equal(y,dV);
    const int m=(int)y.size(), h=(int)y[0].size();
    long long S=0;
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j)
        S += (long long)std::llabs((long long)y[i][j+1]-y[i][j]) * W.wH[i][j];
    for(int i=0;i+1<m;++i) for(int j=0;j<h;++j)
        S += (long long)std::llabs((long long)y[i+1][j]-y[i][j]) * W.wV[i][j];
    return S * dV;
}
bool lightOCShortest::check_OC_lin(const vector<vector<int>>& y){
    const int m=(int)y.size(), h=(int)y[0].size();
    for(int i=0;i<m;++i) for(int j=0;j+1<h;++j) if(y[i][j] > y[i][j+1]) return false;
    for(int j=0;j<h;++j) for(int i=0;i+1<m;++i) if(y[i][j] > y[i+1][j]) return false;
    return true;
}

// ===== ceil log2
static inline int ceil_log2_int(int x){ int r=0, v=x-1; while(v>0){ v>>=1; ++r; } return std::max(1,r); }

// ===== phi build:  phi = (up + left) - (down + right), then * dV =====
void lightOCShortest::build_phi(int m,int h){
    phi_.assign(m, std::vector<long long>(h, 0));
    auto H_at=[&](int I,int J)->long long{
        if(!cfg_.W.enabled) return (0<=I && I<m && 0<=J && J<h-1)? 1LL:0LL;
        return (0<=I && I<m && 0<=J && J<h-1)? cfg_.W.wH[I][J] : 0LL;
    };
    auto V_at=[&](int I,int J)->long long{
        if(!cfg_.W.enabled) return (0<=I && I<m-1 && 0<=J && J<h)? 1LL:0LL;
        return (0<=I && I<m-1 && 0<=J && J<h)? cfg_.W.wV[I][J] : 0LL;
    };

    for(int i=0;i<m;++i) for(int j=0;j<h;++j){
        long long up   = V_at(i-1,j);
        long long left = H_at(i, j-1);
        long long down = V_at(i, j);
        long long right= H_at(i, j);
        long long v = (up + left) - (down + right);
        phi_[i][j] = v * cfg_.dV;
    }
}

// ===== bit-pack codec =====
void lightOCShortest::init_bitpack_codec(int m,int h){
    seg_ = ceil_log2_int(h+1);
    total_bits_ = m * seg_;
    if     (total_bits_ <=  64) nwords_ = 1;
    else if(total_bits_ <= 128) nwords_ = 2;
    else if(total_bits_ <= 192) nwords_ = 3;
    else if(total_bits_ <= 256) nwords_ = 4;
    else throw std::runtime_error("bit-pack needs >256 bits; reduce m,h or add big-int pack");

    mask_ = (seg_==64)? ~0ull : ((1ull<<seg_)-1ull);
    for(int i=0;i<m;++i) shift_[i] = (m-1-i)*seg_;
}
int lightOCShortest::get_digit(const KeyX& k, int idx) const{
    int sh = shift_[idx], w = sh>>6, off=sh&63;
    if(seg_==64) return (int)k.w[w];
    if(off + seg_ <= 64){
        return (int)((k.w[w] >> off) & mask_);
    }else{
        int low = 64 - off;
        uint64_t lo = (k.w[w] >> off);
        uint64_t hi = (k.w[w+1] & ((1ull<<(seg_-low))-1ull));
        return (int)(lo | (hi << low));
    }
}
void lightOCShortest::set_digit(KeyX& k, int idx, int v) const{
    uint64_t vv=(uint64_t)v; int sh=shift_[idx], w=sh>>6, off=sh&63;
    if(seg_==64){ k.w[w]=vv; return; }
    if(off + seg_ <= 64){
        uint64_t m = mask_ << off;
        k.w[w] = (k.w[w] & ~m) | ((vv & mask_) << off);
    }else{
        int low = 64 - off;
        uint64_t mlo = ((1ull<<low)-1ull) << off;
        uint64_t mhi = (1ull<<(seg_-low))-1ull;
        k.w[w]   = (k.w[w]   & ~mlo) | ((vv & ((1ull<<low)-1ull)) << off);
        k.w[w+1] = (k.w[w+1] & ~mhi) | ( (vv >> low) & mhi );
    }
}
lightOCShortest::KeyX lightOCShortest::inc_digit(const KeyX& k, int idx) const{
    KeyX r=k; int v=get_digit(r,idx); set_digit(r,idx,v+1); return r;
}
lightOCShortest::KeyX lightOCShortest::dec_digit(const KeyX& k, int idx) const{
    KeyX r=k; int v=get_digit(r,idx); set_digit(r,idx,v-1); return r;
}

// ===== forward/backward one side (64-bit) =====
struct KV64 { uint64_t key; long long dist; };

static void forward_to_level64(
    const lightOCShortest& S, int m,int h, uint64_t key0,int l0,int steps,
    const vector<uint64_t>& powB, uint64_t B,
    unordered_map<uint64_t,long long>& out,
    bool progress, int omp_threads)
{
    vector<KV64> cur; cur.reserve(4096);
    cur.push_back({key0,0});
    for(int lvl=0; lvl<steps; ++lvl){
        size_t est = std::max<size_t>(cur.size()*4, 16);
        vector<KV64> next; next.reserve(est);
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(est);

#ifdef _OPENMP
        int TH = omp_threads>0? omp_threads : omp_get_max_threads();
        vector<vector<KV64>> loc(TH);
        vector<unordered_map<uint64_t,int>> lid(TH);
        for(int t=0;t<TH;++t){ loc[t].reserve(est/TH+16); lid[t].reserve(est/TH+16); }

#pragma omp parallel num_threads(TH)
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            auto &L = loc[tid]; auto &ID = lid[tid];

#pragma omp for schedule(static)
            for(int u=0; u<(int)cur.size(); ++u){
                auto ku=cur[u].key; auto du=cur[u].dist;
                int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCShortest::digit_at64(ku,i,powB,B);
                for(int i=0;i<m;++i){
                    int aii=ai[i]; if(aii>=h) continue;
                    if(i>0 && aii+1>ai[i-1]) continue;
                    auto kv = lightOCShortest::encode_digit_inc64(ku,i,powB);
                    long long c = (long long)(l0 + lvl + 1) * S.phi_at(i, aii);
                    long long nd = du + c;
                    auto it = ID.find(kv);
                    if(it==ID.end()){ int v=(int)L.size(); L.push_back({kv,nd}); ID.emplace(kv,v); }
                    else{ int v=it->second; if(nd<L[v].dist) L[v].dist=nd; }
                }
            }
        }
        for(auto &vec:loc){
            for(auto &kv: vec){
                auto it=nxt_id.find(kv.key);
                if(it==nxt_id.end()){ int v=(int)next.size(); next.push_back(kv); nxt_id.emplace(kv.key,v); }
                else{ int v=it->second; if(kv.dist<next[v].dist) next[v].dist=kv.dist; }
            }
        }
#else
        for(int u=0; u<(int)cur.size(); ++u){
            auto ku=cur[u].key; auto du=cur[u].dist;
            int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCShortest::digit_at64(ku,i,powB,B);
            for(int i=0;i<m;++i){
                int aii=ai[i]; if(aii>=h) continue;
                if(i>0 && aii+1>ai[i-1]) continue;
                auto kv = lightOCShortest::encode_digit_inc64(ku,i,powB);
                long long c = (long long)(l0 + lvl + 1) * S.phi_at(i, aii);
                long long nd = du + c;
                auto it = nxt_id.find(kv);
                if(it==nxt_id.end()){ int v=(int)next.size(); next.push_back({kv,nd}); nxt_id.emplace(kv, v); }
                else{ int v=it->second; if(nd<next[v].dist) next[v].dist=nd; }
            }
        }
#endif
        if(progress && ((l0+lvl)%16==0))
            cout << "[DP-fwd] L="<<(l0+lvl)<<" states="<<cur.size()<<" -> "<<next.size()<<endl;
        cur.swap(next);
    }
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv:cur) out.emplace(kv.key, kv.dist);
}

static void backward_to_level64(
    const lightOCShortest& S, int m,int h, uint64_t keyN,int lN,int steps,
    const vector<uint64_t>& powB, uint64_t B,
    unordered_map<uint64_t,long long>& out,
    bool progress, int omp_threads)
{
    vector<KV64> cur; cur.reserve(4096);
    cur.push_back({keyN,0});
    for(int t=0; t<steps; ++t){
        int Lcur = lN - t;
        size_t est = std::max<size_t>(cur.size()*4, 16);
        vector<KV64> next; next.reserve(est);
        unordered_map<uint64_t,int> nxt_id; nxt_id.reserve(est);

#ifdef _OPENMP
        int TH = omp_threads>0? omp_threads : omp_get_max_threads();
        vector<vector<KV64>> loc(TH);
        vector<unordered_map<uint64_t,int>> lid(TH);
        for(int k=0;k<TH;++k){ loc[k].reserve(est/TH+16); lid[k].reserve(est/TH+16); }

#pragma omp parallel num_threads(TH)
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            auto &L = loc[tid]; auto &ID = lid[tid];

#pragma omp for schedule(static)
            for(int u=0; u<(int)cur.size(); ++u){
                auto ku=cur[u].key; auto du=cur[u].dist;
                int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCShortest::digit_at64(ku,i,powB,B);
                for(int i=0;i<m;++i){
                    int aii=ai[i]; if(aii<=0) continue;
                    if(i<m-1 && (aii-1)<ai[i+1]) continue;
                    auto kv = lightOCShortest::encode_digit_dec64(ku,i,powB);
                    int jprev = aii - 1;
                    long long c = (long long)Lcur * S.phi_at(i, jprev);
                    long long nd = du + c;
                    auto it = ID.find(kv);
                    if(it==ID.end()){ int v=(int)L.size(); L.push_back({kv,nd}); ID.emplace(kv,v); }
                    else{ int v=it->second; if(nd<L[v].dist) L[v].dist=nd; }
                }
            }
        }
        for(auto &vec:loc){
            for(auto &kv: vec){
                auto it=nxt_id.find(kv.key);
                if(it==nxt_id.end()){ int v=(int)next.size(); next.push_back(kv); nxt_id.emplace(kv.key,v); }
                else{ int v=it->second; if(kv.dist<next[v].dist) next[v].dist=kv.dist; }
            }
        }
#else
        for(int u=0; u<(int)cur.size(); ++u){
            auto ku=cur[u].key; auto du=cur[u].dist;
            int ai[32]; for(int i=0;i<m;++i) ai[i]=lightOCShortest::digit_at64(ku,i,powB,B);
            for(int i=0;i<m;++i){
                int aii=ai[i]; if(aii<=0) continue;
                if(i<m-1 && (aii-1)<ai[i+1]) continue;
                auto kv = lightOCShortest::encode_digit_dec64(ku,i,powB);
                int jprev = aii - 1;
                long long c = (long long)Lcur * S.phi_at(i, jprev);
                long long nd = du + c;
                auto it = nxt_id.find(kv);
                if(it==nxt_id.end()){ int v=(int)next.size(); next.push_back({kv,nd}); nxt_id.emplace(kv, v); }
                else{ int v=it->second; if(nd<next[v].dist) next[v].dist=nd; }
            }
        }
#endif
        if(progress && ((lN-t)%16==0))
            cout << "[DP-bwd] L="<<(lN-t)<<" states="<<cur.size()<<" -> "<<next.size()<<endl;
        cur.swap(next);
    }
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv:cur) out.emplace(kv.key, kv.dist);
}

// ===== forward/backward (KeyX) =====
struct KVX  { lightOCShortest::KeyX key; long long dist; };

static void forward_to_levelX(
    const lightOCShortest& S, int m,int h, lightOCShortest::KeyX key0,int l0,int steps,
    std::unordered_map<lightOCShortest::KeyX,long long,lightOCShortest::KeyXHasher>& out,
    bool progress, int omp_threads)
{
    using KeyX = lightOCShortest::KeyX;
    using Hasher = lightOCShortest::KeyXHasher;
    std::vector<KVX> cur; cur.reserve(4096);
    cur.push_back({key0,0});
    for(int lvl=0; lvl<steps; ++lvl){
        size_t est = std::max<size_t>(cur.size()*4, 16);
        std::vector<KVX> next; next.reserve(est);
        std::unordered_map<KeyX,int,Hasher> nxt_id; nxt_id.reserve(est);

#ifdef _OPENMP
        int TH = omp_threads>0? omp_threads : omp_get_max_threads();
        std::vector<std::vector<KVX>> loc(TH);
        std::vector<std::unordered_map<KeyX,int,Hasher>> lid(TH);
        for(int t=0;t<TH;++t){ loc[t].reserve(est/TH+16); lid[t].reserve(est/TH+16); }

#pragma omp parallel num_threads(TH)
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            auto &L = loc[tid]; auto &ID = lid[tid];

#pragma omp for schedule(static)
            for(int u=0; u<(int)cur.size(); ++u){
                auto ku=cur[u].key; auto du=cur[u].dist;
                int ai[32]; for(int i=0;i<m;++i) ai[i]=S.get_digit(ku,i);
                for(int i=0;i<m;++i){
                    int aii=ai[i]; if(aii>=h) continue;
                    if(i>0 && aii+1>ai[i-1]) continue;
                    auto kv = S.inc_digit(ku,i);
                    long long c = (long long)(l0 + lvl + 1) * S.phi_at(i, aii);
                    long long nd = du + c;
                    auto it = ID.find(kv);
                    if(it==ID.end()){ int v=(int)L.size(); L.push_back({kv,nd}); ID.emplace(kv,v); }
                    else{ int v=it->second; if(nd<L[v].dist) L[v].dist=nd; }
                }
            }
        }
        for(auto &vec:loc){
            for(auto &kv: vec){
                auto it=nxt_id.find(kv.key);
                if(it==nxt_id.end()){ int v=(int)next.size(); next.push_back(kv); nxt_id.emplace(kv.key,v); }
                else{ int v=it->second; if(kv.dist<next[v].dist) next[v].dist=kv.dist; }
            }
        }
#else
        for(int u=0; u<(int)cur.size(); ++u){
            auto ku=cur[u].key; auto du=cur[u].dist;
            int ai[32]; for(int i=0;i<m;++i) ai[i]=S.get_digit(ku,i);
            for(int i=0;i<m;++i){
                int aii=ai[i]; if(aii>=h) continue;
                if(i>0 && aii+1>ai[i-1]) continue;
                auto kv = S.inc_digit(ku,i);
                long long c = (long long)(l0 + lvl + 1) * S.phi_at(i, aii);
                long long nd = du + c;
                auto it = nxt_id.find(kv);
                if(it==nxt_id.end()){ int v=(int)next.size(); next.push_back({kv,nd}); nxt_id.emplace(kv, v); }
                else{ int v=it->second; if(nd<next[v].dist) next[v].dist=nd; }
            }
        }
#endif
        if(progress && ((l0+lvl)%16==0))
            std::cout << "[DP-fwdX] L="<<(l0+lvl)<<" states="<<cur.size()<<" -> "<<next.size()<<std::endl;
        cur.swap(next);
    }
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv:cur) out.emplace(kv.key, kv.dist);
}

static void backward_to_levelX(
    const lightOCShortest& S, int m,int h, lightOCShortest::KeyX keyN,int lN,int steps,
    std::unordered_map<lightOCShortest::KeyX,long long,lightOCShortest::KeyXHasher>& out,
    bool progress, int omp_threads)
{
    using KeyX = lightOCShortest::KeyX;
    using Hasher = lightOCShortest::KeyXHasher;
    std::vector<KVX> cur; cur.reserve(4096);
    cur.push_back({keyN,0});
    for(int t=0; t<steps; ++t){
        int Lcur = lN - t;
        size_t est = std::max<size_t>(cur.size()*4, 16);
        std::vector<KVX> next; next.reserve(est);
        std::unordered_map<KeyX,int,Hasher> nxt_id; nxt_id.reserve(est);

#ifdef _OPENMP
        int TH = omp_threads>0? omp_threads : omp_get_max_threads();
        std::vector<std::vector<KVX>> loc(TH);
        std::vector<std::unordered_map<KeyX,int,Hasher>> lid(TH);
        for(int k=0;k<TH;++k){ loc[k].reserve(est/TH+16); lid[k].reserve(est/TH+16); }

#pragma omp parallel num_threads(TH)
        {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            auto &L = loc[tid]; auto &ID = lid[tid];

#pragma omp for schedule(static)
            for(int u=0; u<(int)cur.size(); ++u){
                auto ku=cur[u].key; auto du=cur[u].dist;
                int ai[32]; for(int i=0;i<m;++i) ai[i]=S.get_digit(ku,i);
                for(int i=0;i<m;++i){
                    int aii=ai[i]; if(aii<=0) continue;
                    if(i<m-1 && (aii-1)<ai[i+1]) continue;
                    auto kv = S.dec_digit(ku,i);
                    int jprev = aii - 1;
                    long long c = (long long)Lcur * S.phi_at(i, jprev);
                    long long nd = du + c;
                    auto it = ID.find(kv);
                    if(it==ID.end()){ int v=(int)L.size(); L.push_back({kv,nd}); ID.emplace(kv,v); }
                    else{ int v=it->second; if(nd<L[v].dist) L[v].dist=nd; }
                }
            }
        }
        for(auto &vec:loc){
            for(auto &kv: vec){
                auto it=nxt_id.find(kv.key);
                if(it==nxt_id.end()){ int v=(int)next.size(); next.push_back(kv); nxt_id.emplace(kv.key,v); }
                else{ int v=it->second; if(kv.dist<next[v].dist) next[v].dist=kv.dist; }
            }
        }
#else
        for(int u=0; u<(int)cur.size(); ++u){
            auto ku=cur[u].key; auto du=cur[u].dist;
            int ai[32]; for(int i=0;i<m;++i) ai[i]=S.get_digit(ku,i);
            for(int i=0;i<m;++i){
                int aii=ai[i]; if(aii<=0) continue;
                if(i<m-1 && (aii-1)<ai[i+1]) continue;
                auto kv = S.dec_digit(ku,i);
                int jprev = aii - 1;
                long long c = (long long)Lcur * S.phi_at(i, jprev);
                long long nd = du + c;
                auto it = nxt_id.find(kv);
                if(it==nxt_id.end()){ int v=(int)next.size(); next.push_back({kv,nd}); nxt_id.emplace(kv, v); }
                else{ int v=it->second; if(nd<next[v].dist) next[v].dist=nd; }
            }
        }
#endif
        if(progress && ((lN-t)%16==0))
            std::cout << "[DP-bwdX] L="<<(lN-t)<<" states="<<cur.size()<<" -> "<<next.size()<<std::endl;
        cur.swap(next);
    }
    out.clear(); out.reserve(cur.size()*2+16);
    for(auto &kv:cur) out.emplace(kv.key, kv.dist);
}

// ===== reconstruct (Hirschberg) =====
void lightOCShortest::reconstruct_hirschberg64(
    int m,int h,long long dV,
    const std::vector<uint64_t>& powB, uint64_t B, bool progress, int omp_threads,
    int l_lo, uint64_t key_lo, int l_hi, uint64_t key_hi,
    vector<vector<int>>& y)
{
    int len = l_hi - l_lo;
    if(len==0) return;
    if(len==1){
        int sel_i=-1;
        for(int i=0;i<m;++i){
            int a_lo = digit_at64(key_lo,i,powB,B);
            int a_hi = digit_at64(key_hi,i,powB,B);
            if(a_hi == a_lo + 1){ sel_i=i; break; }
        }
        if(sel_i<0) throw std::runtime_error("Base step mismatch (64)");
        int j0 = digit_at64(key_lo, sel_i, powB, B);
        y[sel_i][j0] = l_lo + 1;
        return;
    }
    int l_mid = l_lo + len/2;

    std::unordered_map<uint64_t,long long> F, Bmap;
    forward_to_level64(*this, m,h, key_lo, l_lo, l_mid-l_lo, powB, B, F, progress, omp_threads);
    backward_to_level64(*this,m,h, key_hi, l_hi, l_hi-l_mid, powB, B, Bmap, progress, omp_threads);

    uint64_t key_mid=0; long long best=std::numeric_limits<long long>::max();
    if(F.size()<=Bmap.size()){
        for(auto &kv:F){
            auto it=Bmap.find(kv.first);
            if(it==Bmap.end()) continue;
            long long v=kv.second + it->second;
            if(v<best){ best=v; key_mid=kv.first; }
        }
    }else{
        for(auto &kv:Bmap){
            auto it=F.find(kv.first);
            if(it==F.end()) continue;
            long long v=kv.second + it->second;
            if(v<best){ best=v; key_mid=kv.first; }
        }
    }
    if(best==std::numeric_limits<long long>::max())
        throw std::runtime_error("No midpoint on shortest path (64)");

    reconstruct_hirschberg64(m,h,dV,powB,B,progress,omp_threads, l_lo,key_lo,l_mid,key_mid,y);
    reconstruct_hirschberg64(m,h,dV,powB,B,progress,omp_threads, l_mid,key_mid,l_hi,key_hi,y);
}

void lightOCShortest::reconstruct_hirschbergX(
    int m,int h,long long dV, bool progress, int omp_threads,
    int l_lo, KeyX key_lo, int l_hi, KeyX key_hi,
    vector<vector<int>>& y)
{
    int len = l_hi - l_lo;
    if(len==0) return;
    if(len==1){
        int sel_i=-1;
        for(int i=0;i<m;++i){
            int a_lo = get_digit(key_lo,i);
            int a_hi = get_digit(key_hi,i);
            if(a_hi == a_lo + 1){ sel_i=i; break; }
        }
        if(sel_i<0) throw std::runtime_error("Base step mismatch (KeyX)");
        int j0 = get_digit(key_lo, sel_i);
        y[sel_i][j0] = l_lo + 1;
        return;
    }
    int l_mid = l_lo + len/2;

    std::unordered_map<KeyX,long long,KeyXHasher> F, Bmap;
    forward_to_levelX(*this, m,h, key_lo, l_lo, l_mid-l_lo, F, progress, omp_threads);
    backward_to_levelX(*this,m,h, key_hi, l_hi, l_hi-l_mid, Bmap, progress, omp_threads);

    KeyX key_mid{}; long long best=std::numeric_limits<long long>::max();
    if(F.size()<=Bmap.size()){
        for(auto &kv:F){
            auto it=Bmap.find(kv.first);
            if(it==Bmap.end()) continue;
            long long v=kv.second + it->second;
            if(v<best){ best=v; key_mid=kv.first; }
        }
    }else{
        for(auto &kv:Bmap){
            auto it=F.find(kv.first);
            if(it==F.end()) continue;
            long long v=kv.second + it->second;
            if(v<best){ best=v; key_mid=kv.first; }
        }
    }
    if(best==std::numeric_limits<long long>::max())
        throw std::runtime_error("No midpoint on shortest path (KeyX)");

    reconstruct_hirschbergX(m,h,dV,progress,omp_threads, l_lo,key_lo,l_mid,key_mid,y);
    reconstruct_hirschbergX(m,h,dV,progress,omp_threads, l_mid,key_mid,l_hi,key_hi,y);
}

// ===== solve =====
OCShortestResult lightOCShortest::solve(int m, int h){
    const int n_full = m*h;
    const int n = (cfg_.max_steps>0 && cfg_.max_steps<n_full) ? cfg_.max_steps : n_full;

    build_phi(m,h);

    std::vector<std::vector<int>> y(m, std::vector<int>(h, 0));

    if(!cfg_.use_pow2_pack){
        // small: 64-bit base-(h+1)
        const uint64_t B = (uint64_t)h + 1u;
        vector<uint64_t> powB(m);
        powB[m-1]=1ull;
        for(int k=m-2;k>=0;--k){
            __uint128_t tmp=(__uint128_t)powB[k+1]*(__uint128_t)B;
            if(tmp > (__uint128_t)std::numeric_limits<uint64_t>::max())
                throw std::runtime_error("64-bit key overflow; run with --pow2=1");
            powB[k]=(uint64_t)tmp;
        }
        uint64_t key0=0ull, keyN=0ull;
        for(int i=0;i<m;++i) keyN += powB[i] * (uint64_t)h;

        if(cfg_.verbose){
            cout << "[OC-Shortest/64] m="<<m<<" h="<<h<<" levels="<<n;
#ifdef _OPENMP
            cout << " threads=" << (cfg_.omp_threads>0?cfg_.omp_threads:omp_get_max_threads());
#endif
            cout << endl;
        }
        reconstruct_hirschberg64(m,h,cfg_.dV,powB,B,cfg_.progress,cfg_.omp_threads, 0,key0, n,keyN, y);
    }else{
        // big: 2^k bit-pack
        init_bitpack_codec(m,h);
        KeyX key0{}, keyN{};
        for(int i=0;i<m;++i) set_digit(keyN,i,h);

        if(cfg_.verbose){
            cout << "[OC-Shortest/pow2] m="<<m<<" h="<<h
                 << " SEG="<<seg_<<" bits="<<total_bits_<<" words="<<nwords_;
#ifdef _OPENMP
            cout << " threads=" << (cfg_.omp_threads>0?cfg_.omp_threads:omp_get_max_threads());
#endif
            cout << endl;
        }
        reconstruct_hirschbergX(m,h,cfg_.dV,cfg_.progress,cfg_.omp_threads, 0,key0, n,keyN, y);
    }

    // report HPWL with dV=1 (phi already carries dV)
    long long total = cfg_.W.enabled ? hpwl_sum_weighted(y, cfg_.W, 1)
                                     : hpwl_sum_equal(y, 1);
    OCShortestResult R{m,h,n,cfg_.dV, std::move(y), total};
    return R;
}
