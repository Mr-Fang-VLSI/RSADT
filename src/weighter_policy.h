#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

template<typename T> using Mat = std::vector<std::vector<T>>;

struct WeighterReport {
    long long WNS = 0;         // min(T - L_e) under the T passed to update
    long long TNS = 0;         // sum of violations
    int picked_edges = 0;      // edges updated this round (policy-defined)
};

struct WeighterParams {
    // Export/clip
    double export_scale = 1000.0;
    double cap_max = 32.0;

    // Momentum (adapter)
    double alpha = 0.5;
    double eta   = 0.5;
    double top_ratio = 0.10;

    // Lagrangian
    double rho = 0.3;

    // IRL1
    double lam = 2.0;
    double delta = 2.0;

    // Soft-threshold
    double gamma = 0.2, beta = 2.0;
};

class IWeighterPolicy {
public:
    virtual ~IWeighterPolicy() = default;
    virtual std::string name() const = 0;
    virtual void set_initial(const Mat<long long>& wH, const Mat<long long>& wV, long long scale) = 0;
    virtual void export_ll(Mat<long long>& wH_out, Mat<long long>& wV_out) = 0;
    virtual WeighterReport update_from_layout(const Mat<int>& y, long long T, bool pick_by_ratio=true) = 0;
    virtual void freeze_if_feasible(bool enable) = 0;

    // 新增：清空“惯性/动量”缓冲，但保留当前权重（默认 no-op）
    virtual void reset_inertia() {}
};

// factories
IWeighterPolicy* MakeMomentumPolicy(int m, int h, const WeighterParams& P);
IWeighterPolicy* MakeLagrangePolicy(int m, int h, const WeighterParams& P);
IWeighterPolicy* MakeIRL1Policy(int m, int h, const WeighterParams& P);
IWeighterPolicy* MakeSoftPolicy(int m, int h, const WeighterParams& P);
IWeighterPolicy* MakeFreezeDecorator(IWeighterPolicy* inner);
