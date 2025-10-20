// rsad_cp.cpp
// Baseline CP-SAT model for chain placement with OC + Max-T + AllDifferent
// Objective: boundary-decomposed linear form (R-SAD Lemma 1, Eq.(2))
// Author: (c) 2025

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <cassert>

#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;

using operations_research::sat::CpModelBuilder;
using operations_research::sat::CpSolverResponse;
using operations_research::sat::CpSolverStatus;
using operations_research::sat::CpSolverStatus_Name;
using operations_research::sat::Domain;
using operations_research::sat::IntVar;
using operations_research::sat::LinearExpr;
using operations_research::sat::Model;
using operations_research::sat::SatParameters;
using operations_research::sat::Solve;
using operations_research::sat::SolutionIntegerValue;

struct SolveResult {
    CpSolverStatus status;
    int64_t objective_value = 0;
    vector<vector<int64_t>> y;  // m x h assignment
    int64_t vertical_hpwl = 0;  // computed from solution via boundary formula
};

class GridChainPlacerCP {
public:
    GridChainPlacerCP(int m, int h, int T)
        : m_(m), h_(h), n_(m * h), T_(T) {
        if (m_ <= 0 || h_ <= 0) throw std::invalid_argument("m,h must be positive.");
        if (T_ < 1) throw std::invalid_argument("T must be >= 1.");
        weights_ = BuildBoundaryWeights(/*drop_corner_constants=*/true);
    }

    // 可替换为自定义权重（要求与y同shape）
    void SetWeights(const vector<vector<int>>& w) {
        if ((int)w.size() != m_) throw std::invalid_argument("weights rows != m");
        for (int i = 0; i < m_; ++i) if ((int)w[i].size() != h_) throw std::invalid_argument("weights cols != h");
        weights_ = w;
    }

    // 运行求解（time_limit_sec=0 则不设时限）
    SolveResult SolveOnce(double time_limit_sec = 30.0, bool verbose_log = false, bool use_pairwise_all_diff = true) const {
        CpModelBuilder cp;
        // 变量 y[i][j] ∈ [LB, UB]，LB/UB来自OC的最紧链式界
        vector<vector<IntVar>> Y(m_, vector<IntVar>(h_));
        for (int i = 0; i < m_; ++i) {
            for (int j = 0; j < h_; ++j) {
                int64_t lb = 1 + i + j; // 1-based最小层级
                int64_t ub = n_ - ((m_ - 1 - i) + (h_ - 1 - j));
                Y[i][j] = cp.NewIntVar(Domain(lb, ub), "y_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        // OC: strict increasing along rows/cols
        for (int i = 0; i < m_; ++i) {
            for (int j = 0; j < h_; ++j) {
                if (i + 1 < m_) {
                    cp.AddLessOrEqual(LinearExpr(Y[i][j]) + 1, Y[i + 1][j]);     // y[i+1][j] >= y[i][j] + 1
                    cp.AddLessOrEqual(Y[i + 1][j], LinearExpr(Y[i][j]) + T_);    // Max-T
                }
                if (j + 1 < h_) {
                    cp.AddLessOrEqual(LinearExpr(Y[i][j]) + 1, Y[i][j + 1]);     // y[i][j+1] >= y[i][j] + 1
                    cp.AddLessOrEqual(Y[i][j + 1], LinearExpr(Y[i][j]) + T_);    // Max-T
                }
            }
        }

        // 占位唯一性：AllDifferent over all cells.
        // CP-SAT C++ 也支持 AddAllDifferent，但为兼容不同版本，这里用成对!=（O(n^2)）
        vector<IntVar> flat; flat.reserve(n_);
        for (int i = 0; i < m_; ++i) for (int j = 0; j < h_; ++j) flat.push_back(Y[i][j]);
        if (use_pairwise_all_diff) {
            for (int a = 0; a < n_; ++a) {
                for (int b = a + 1; b < n_; ++b) {
                    cp.AddNotEqual(flat[a], flat[b]);
                }
            }
        } else {
            cp.AddAllDifferent(flat); // 若你的 OR-Tools 版本支持，亦可使用此句
        }

        // 可选：显式锚点（在OC+AllDifferent下本身就会成立，但能增强传播）
        cp.AddEquality(Y[0][0], 1);
        cp.AddEquality(Y[m_ - 1][h_ - 1], n_);

        // 目标：∑ w[i][j] * y[i][j]   （R-SAD Lemma 1 式(2)的线性化）
        LinearExpr obj;
        for (int i = 0; i < m_; ++i)
            for (int j = 0; j < h_; ++j)
                if (weights_[i][j] != 0) obj += weights_[i][j] * Y[i][j];
        cp.Minimize(obj);

        // 参数
        Model mdl;
        SatParameters p;
        if (time_limit_sec > 0.0) p.set_max_time_in_seconds(time_limit_sec);
        p.set_num_search_workers(8);                 // 并行（按需调整）
        p.set_log_search_progress(verbose_log);      // 控制日志
        mdl.Add(operations_research::sat::NewSatParameters(p));

        // 求解
        const CpSolverResponse resp = Solve(cp.Build(), &mdl);

        SolveResult res;
        res.status = resp.status();
        if (resp.status() == CpSolverStatus::OPTIMAL || resp.status() == CpSolverStatus::FEASIBLE) {
            res.y.assign(m_, vector<int64_t>(h_, -1));
            for (int i = 0; i < m_; ++i)
                for (int j = 0; j < h_; ++j)
                    res.y[i][j] = SolutionIntegerValue(resp, Y[i][j]);
            res.objective_value = resp.objective_value();

            // 计算“纵向HPWL边界和”（与目标同形，用于 sanity check）
            int64_t hpwl_v = 0;
            for (int j = 0; j < h_; ++j) hpwl_v += res.y[m_-1][j] - res.y[0][j];
            for (int i = 0; i < m_; ++i) hpwl_v += res.y[i][h_-1] - res.y[i][0];
            res.vertical_hpwl = hpwl_v;
        }
        return res;
    }

    // 生成边界权重（-1/0/+1），默认去掉常数项（(1,1)与(m,h)）
    vector<vector<int>> BuildBoundaryWeights(bool drop_corner_constants = true) const {
        vector<vector<int>> w(m_, vector<int>(h_, 0));
        for (int i = 0; i < m_; ++i) {
            for (int j = 0; j < h_; ++j) {
                int val = 0;
                if (i == m_ - 1) val += 1;  // + y(m, j)
                if (i == 0)      val -= 1;  // - y(1, j)
                if (j == h_ - 1) val += 1;  // + y(i, h)
                if (j == 0)      val -= 1;  // - y(i, 1)
                w[i][j] = val;
            }
        }
        if (drop_corner_constants) {
            // (1,1)与(m,h)的 ±2 仅是常数项（在OC下 y(1,1)=1, y(m,h)=n），可安全置零
            w[0][0] = 0;
            w[m_-1][h_-1] = 0;
        }
        return w;
    }

private:
    int m_, h_, n_, T_;
    vector<vector<int>> weights_;
};

// --------------------- demo main ---------------------
static void PrintSolution(const SolveResult& r) {
    cout << "Status: " << CpSolverStatus_Name(r.status) << "\n";
    if (r.status != operations_research::sat::CpSolverStatus::OPTIMAL &&
        r.status != operations_research::sat::CpSolverStatus::FEASIBLE) return;

    cout << "Objective (boundary linear form) = " << r.objective_value << "\n";
    cout << "Vertical HPWL (computed)        = " << r.vertical_hpwl << "\n";
    cout << "Assignment y[i][j] (1..n), row-major view:\n";
    for (const auto& row : r.y) {
        for (auto v : row) cout << v << "\t";
        cout << "\n";
    }
}

int main(int argc, char** argv) {
    // 示例：m h T 可从命令行传入；默认 8x8, T=8（接近 min(m,h)，几乎退化为仅OC）
    int m = 8, h = 8, T = 8;
    if (argc >= 3) { m = std::stoi(argv[1]); h = std::stoi(argv[2]); }
    if (argc >= 4) { T = std::stoi(argv[3]); }
    cout << "Solving m=" << m << ", h=" << h << ", T=" << T << " ...\n";

    try {
        GridChainPlacerCP solver(m, h, T);
        // 也可以自定义权重（例如只优化列边界、或附加其他线性项）
        // auto w = solver.BuildBoundaryWeights(true); solver.SetWeights(w);

        auto res = solver.SolveOnce(/*time_limit_sec=*/30.0, /*verbose_log=*/false, /*use_pairwise_all_diff=*/true);
        PrintSolution(res);
        if (res.status != operations_research::sat::CpSolverStatus::OPTIMAL &&
            res.status != operations_research::sat::CpSolverStatus::FEASIBLE) {
            cerr << "No feasible solution under given T. Try increasing T.\n";
            return 2;
        }
    } catch (const std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
