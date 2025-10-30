#include "weighter_policy.h"
#include "momentum_weighter.h"
#include <numeric>
#include <cassert>
#include <cmath>

// ---------- utils ----------
template<typename T> using MatT = std::vector<std::vector<T>>;
static inline long long llround_clip(double x){ return (long long) std::llround(x); }

static void clip_and_scale(MatT<double>& W, double scale, double cap){
    const double lo = scale, hi = cap*scale;
    for(auto& r: W) for(auto& v: r){ v = std::min(std::max(v*scale, lo), hi); }
}

static void to_ll(const MatT<double>& Wd, MatT<long long>& Wll){
    Wll.assign(Wd.size(), {});
    for(size_t i=0;i<Wd.size();++i){
        Wll[i].resize(Wd[i].size());
        for(size_t j=0;j<Wd[i].size();++j) Wll[i][j] = llround_clip(Wd[i][j]);
    }
}

struct ScanOut {
    WeighterReport rep;
    MatT<long long> LH, LV;
};

static ScanOut scan_edges(const Mat<int>& y, int m, int h, long long T){
    ScanOut S; S.LH.assign(m, std::vector<long long>(h-1,0));
    S.LV.assign(m-1, std::vector<long long>(h,0));
    long long WNS = (long long)1e15, TNS = 0;
    for(int i=0;i<m;++i)for(int j=0;j<h-1;++j){
        long long L = std::llabs((long long)y[i][j+1] - (long long)y[i][j]); S.LH[i][j]=L;
        long long s = T - L; WNS = std::min(WNS, s); if(s<0) TNS += -s;
    }
    for(int i=0;i<m-1;++i)for(int j=0;j<h;++j){
        long long L = std::llabs((long long)y[i+1][j] - (long long)y[i][j]); S.LV[i][j]=L;
        long long s = T - L; WNS = std::min(WNS, s); if(s<0) TNS += -s;
    }
    S.rep.WNS = (WNS==(long long)1e15? 0: WNS);
    S.rep.TNS = TNS;
    return S;
}

// ---------- Momentum adapter（与你现有实现一致） ----------
class MomentumPolicy final : public IWeighterPolicy {
    MomentumWeighter W;
    Mat<long long> last_wH, last_wV;
    bool frozen=false;
public:
    MomentumPolicy(int m, int h, const WeighterParams& P) : W(m,h, MomentumWeighter::Params{
        P.alpha, P.eta, P.top_ratio, P.cap_max, P.export_scale
    }) {}
    std::string name() const override { return "momentum"; }
    void set_initial(const Mat<long long>& wH, const Mat<long long>& wV, long long scale) override {
        W.set_initial(wH, wV, scale);
    }
    void export_ll(Mat<long long>& wH_out, Mat<long long>& wV_out) override {
        if (frozen){ wH_out=last_wH; wV_out=last_wV; return; }
        W.export_ll(wH_out, wV_out); last_wH=wH_out; last_wV=wV_out;
    }
    WeighterReport update_from_layout(const Mat<int>& y, long long T, bool pick_by_ratio) override {
        if (frozen){ auto S=scan_edges(y,(int)y.size(),(int)y[0].size(),T); S.rep.picked_edges=0; return S.rep; }
        auto R = W.update_from_layout(y, T, pick_by_ratio);
        return {R.WNS, R.TNS, R.picked_edges};
    }
    void freeze_if_feasible(bool en) override { frozen=en; }
    void reset_inertia() override { /* no-op */ }
};
IWeighterPolicy* MakeMomentumPolicy(int m, int h, const WeighterParams& P){ return new MomentumPolicy(m,h,P); }

// ---------- Lagrangian（与你现有实现一致） ----------
class LagrangePolicy final : public IWeighterPolicy {
    int m,h; WeighterParams P;
    MatT<double> lamH, lamV; MatT<double> wH, wV; bool frozen=false;
public:
    LagrangePolicy(int m_, int h_, const WeighterParams& P_) : m(m_), h(h_), P(P_){
        lamH.assign(m, std::vector<double>(h-1,0.0));
        lamV.assign(m-1, std::vector<double>(h,0.0));
        wH.assign(m, std::vector<double>(h-1,1.0));
        wV.assign(m-1, std::vector<double>(h,1.0));
    }
    std::string name() const override { return "lagrange"; }
    void set_initial(const Mat<long long>&, const Mat<long long>&, long long) override {}
    void export_ll(Mat<long long>& outH, Mat<long long>& outV) override {
        for(int i=0;i<m;++i)for(int j=0;j<h-1;++j) wH[i][j] = 1.0 + lamH[i][j];
        for(int i=0;i<m-1;++i)for(int j=0;j<h;++j)   wV[i][j] = 1.0 + lamV[i][j];
        clip_and_scale(wH, P.export_scale, P.cap_max);
        clip_and_scale(wV, P.export_scale, P.cap_max);
        to_ll(wH,outH); to_ll(wV,outV);
    }
    WeighterReport update_from_layout(const Mat<int>& y, long long T, bool) override {
        auto S=scan_edges(y,m,h,T);
        int picked=0;
        if(!frozen){
            for(int i=0;i<m;++i)for(int j=0;j<h-1;++j){ long long L=S.LH[i][j]; if(L>T){ lamH[i][j]+=P.rho*(L-T); picked++; } }
            for(int i=0;i<m-1;++i)for(int j=0;j<h;++j){ long long L=S.LV[i][j]; if(L>T){ lamV[i][j]+=P.rho*(L-T); picked++; } }
        }
        S.rep.picked_edges=picked; return S.rep;
    }
    void freeze_if_feasible(bool en) override { frozen=en; }
    void reset_inertia() override { /* no-op */ }
};
IWeighterPolicy* MakeLagrangePolicy(int m, int h, const WeighterParams& P){ return new LagrangePolicy(m,h,P); }

// ---------- IRL1⁺（改进版：单调非减 + 阶段式 λ + 动态 δ） ----------
class IRL1Policy final : public IWeighterPolicy {
    int m,h; WeighterParams P;
    MatT<double> wH, wV;   // 当前（未缩放）权重
    bool frozen=false;

    // 阶段控制：当主程序把 T′ 收紧（传入的 T 变小）时，λ_eff 按倍乘放大
    long long last_T_seen = -1;
    double lam_eff;             // 当前阶段的等效 λ
    const double lam_growth = 1.25;        // 每收紧一次 T′，λ 放大倍数
    const double lam_eff_cap_factor =  (/* 上限：不超过导出 cap */ 1.0); // 先用 1.0；下方用 cap 约束最终权重
    const double kappa_delta = 0.25;       // δ_eff = max(δ, κ·T)

public:
    IRL1Policy(int m_, int h_, const WeighterParams& P_) : m(m_), h(h_), P(P_){
        wH.assign(m, std::vector<double>(h-1,1.0));
        wV.assign(m-1, std::vector<double>(h,1.0));
        lam_eff = std::max(1e-9, P.lam);  // 从命令行 lam 启动
    }
    std::string name() const override { return "irl1"; }
    void set_initial(const Mat<long long>&, const Mat<long long>&, long long) override {}

    void export_ll(Mat<long long>& outH, Mat<long long>& outV) override {
        // clip: [scale, cap*scale]
        clip_and_scale(wH, P.export_scale, P.cap_max);
        clip_and_scale(wV, P.export_scale, P.cap_max);
        to_ll(wH,outH); to_ll(wV,outV);
    }

    WeighterReport update_from_layout(const Mat<int>& y, long long T, bool) override {
        auto S=scan_edges(y,m,h,T);
        int picked=0;
        if(!frozen){
            // 1) 阶段判定：若 T（实为 T′）变小，放大 λ_eff
            if (last_T_seen < 0) {
                last_T_seen = T;
            } else if (T < last_T_seen) {
                // 收紧阶段：放大 λ，防止“压到当前 T′ 后不够力再压更紧 T′”
                lam_eff = lam_eff * lam_growth;
                last_T_seen = T;
            }

            // 2) 动态 δ（与 T 同量纲），避免小 T 时过敏
            const double delta_eff = std::max(1.0, std::max(P.delta, kappa_delta * (double)T));

            // 3) 单调非减：仅对违例边更新，w = max(w, 1 + λ_eff * min(1, (L-T)/δ_eff))
            for(int i=0;i<m;++i)for(int j=0;j<h-1;++j){
                const long long L=S.LH[i][j];
                if(L > T){
                    const double v = std::min(1.0, (double)(L - T) / delta_eff);
                    const double target = 1.0 + lam_eff * v;
                    if (target > wH[i][j] + 1e-12) { wH[i][j] = target; picked++; }
                }
            }
            for(int i=0;i<m-1;++i)for(int j=0;j<h;++j){
                const long long L=S.LV[i][j];
                if(L > T){
                    const double v = std::min(1.0, (double)(L - T) / delta_eff);
                    const double target = 1.0 + lam_eff * v;
                    if (target > wV[i][j] + 1e-12) { wV[i][j] = target; picked++; }
                }
            }

            // 注：不对非违例边降权（避免“刚压好又松回去”的反复）
            // λ_eff 的最终作用仍由 export 时 cap 控制到 [scale, cap*scale] 范围内
        }
        S.rep.picked_edges = picked; return S.rep;
    }

    void freeze_if_feasible(bool en) override { frozen=en; }
    void reset_inertia() override { /* IRL1 无动量，保持空实现 */ }
};
IWeighterPolicy* MakeIRL1Policy(int m, int h, const WeighterParams& P){ return new IRL1Policy(m,h,P); }

// ---------- Soft-threshold（你当前的 DreamPlace 风格：log 累积 + α 惯性） ----------
class SoftPolicy final : public IWeighterPolicy {
    int m,h; WeighterParams P;
    MatT<double> logH, logV;     // log w
    MatT<double> dlogH, dlogV;   // 动量(EMA)
    bool frozen=false;

    static inline double phi(double x){ return std::log(1.0 + x); } // smooth
public:
    SoftPolicy(int m_, int h_, const WeighterParams& P_) : m(m_), h(h_), P(P_){
        logH.assign(m, std::vector<double>(h-1, 0.0)); // w=exp(0)=1
        logV.assign(m-1, std::vector<double>(h,   0.0));
        dlogH.assign(m, std::vector<double>(h-1, 0.0));
        dlogV.assign(m-1, std::vector<double>(h,   0.0));
    }
    std::string name() const override { return "soft"; }
    void set_initial(const Mat<long long>&, const Mat<long long>&, long long) override {}

    void export_ll(Mat<long long>& outH, Mat<long long>& outV) override {
        MatT<double> wH(m, std::vector<double>(h-1,1.0));
        MatT<double> wV(m-1, std::vector<double>(h,1.0));
        for(int i=0;i<m;++i)for(int j=0;j<h-1;++j) wH[i][j]=std::exp(logH[i][j]);
        for(int i=0;i<m-1;++i)for(int j=0;j<h;++j)   wV[i][j]=std::exp(logV[i][j]);
        clip_and_scale(wH, P.export_scale, P.cap_max);
        clip_and_scale(wV, P.export_scale, P.cap_max);
        to_ll(wH,outH); to_ll(wV,outV);
    }

    WeighterReport update_from_layout(const Mat<int>& y, long long T, bool) override {
        auto S=scan_edges(y,m,h,T);
        int picked=0;
        if(!frozen){
            const double a = P.alpha;
            const double one_minus_a = 1.0 - a;
            const double eta = P.eta;

            // 水平边
            for(int i=0;i<m;++i)for(int j=0;j<h-1;++j){
                long long L=S.LH[i][j];
                if(L>T){
                    double v = (double)(L-T)/(double)T;         // 违例强度（归一化）
                    double inc = P.gamma * phi(P.beta * v);     // 当前信号
                    dlogH[i][j] = a * dlogH[i][j] + one_minus_a * inc;  // EMA
                    logH[i][j]  += eta * dlogH[i][j];                   // 累积
                    picked++;
                }else{
                    // 不降权：只让动量自然衰减
                    dlogH[i][j] = a * dlogH[i][j];
                }
            }
            // 垂直边
            for(int i=0;i<m-1;++i)for(int j=0;j<h;++j){
                long long L=S.LV[i][j];
                if(L>T){
                    double v = (double)(L-T)/(double)T;
                    double inc = P.gamma * phi(P.beta * v);
                    dlogV[i][j] = a * dlogV[i][j] + one_minus_a * inc;
                    logV[i][j]  += eta * dlogV[i][j];
                    picked++;
                }else{
                    dlogV[i][j] = a * dlogV[i][j];
                }
            }
        }
        S.rep.picked_edges = picked; return S.rep;
    }

    void freeze_if_feasible(bool en) override { frozen=en; }
    void reset_inertia() override {
        for(int i=0;i<m;++i)for(int j=0;j<h-1;++j) dlogH[i][j] = 0.0;
        for(int i=0;i<m-1;++i)for(int j=0;j<h;++j) dlogV[i][j] = 0.0;
    }
};
IWeighterPolicy* MakeSoftPolicy(int m, int h, const WeighterParams& P){ return new SoftPolicy(m,h,P); }

// ---------- Freeze decorator（转发 reset_inertia） ----------
class FreezeDecorator final : public IWeighterPolicy {
    IWeighterPolicy* inner; bool frozen=false; Mat<long long> keepH, keepV;
public:
    explicit FreezeDecorator(IWeighterPolicy* p): inner(p){}
    ~FreezeDecorator(){ delete inner; }
    std::string name() const override { return inner->name()+"+freeze"; }
    void set_initial(const Mat<long long>& wH, const Mat<long long>& wV, long long s) override { inner->set_initial(wH,wV,s); }
    void export_ll(Mat<long long>& wH, Mat<long long>& wV) override {
        if (frozen){ wH=keepH; wV=keepV; return; }
        inner->export_ll(wH,wV); keepH=wH; keepV=wV;
    }
    WeighterReport update_from_layout(const Mat<int>& y, long long T, bool pick_by_ratio) override {
        if (frozen){ auto S=scan_edges(y,(int)y.size(),(int)y[0].size(),T); S.rep.picked_edges=0; return S.rep; }
        return inner->update_from_layout(y,T,pick_by_ratio);
    }
    void freeze_if_feasible(bool en) override { frozen=en; inner->freeze_if_feasible(en); }
    void reset_inertia() override { inner->reset_inertia(); }
};
IWeighterPolicy* MakeFreezeDecorator(IWeighterPolicy* inner){ return new FreezeDecorator(inner); }
