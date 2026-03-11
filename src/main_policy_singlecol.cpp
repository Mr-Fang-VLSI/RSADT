#include "lightOCShortest.h"
#include "weighter_policy.h"
#include "rowswap_ctr_expander.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <cmath>

template<typename T> using Mat = std::vector<std::vector<T>>;
using namespace std;

static bool parse_flag(const string& s, const string& key, string& out){
    auto p=s.find('='); if(p==string::npos) return false;
    if(s.substr(0,p)==key){ out=s.substr(p+1); return true; }
    return false;
}

struct StatTW { long long maxL=0, WNS=0; };
static StatTW compute_maxL_WNS(const Mat<int>& y, int m, int h, long long T){
    StatTW S; long long wns=(long long)1e15; long long maxL=0;
    for(int i=0;i<m;++i)for(int j=0;j<h-1;++j){
        long long L= llabs((long long)y[i][j+1]- (long long)y[i][j]); maxL=max(maxL,L); wns=min(wns, T-L);
    }
    for(int i=0;i<m-1;++i)for(int j=0;j<h;++j){
        long long L= llabs((long long)y[i+1][j]- (long long)y[i][j]); maxL=max(maxL,L); wns=min(wns, T-L);
    }
    S.maxL=maxL; S.WNS=(wns==(long long)1e15? 0: wns); return S;
}

static void write_singlecol_placement(const Mat<int>& y, int m, int h, long long T){
    string fn = "placement"+to_string(m)+"_Col_1_T_"+to_string(T)+".txt";
    ofstream fout(fn);
    fout << "This is a sample " << m << "x" << h << " PE array being placed on 1 DSP column:\n\n";
    for(int i=0;i<m;++i){
        for(int j=0;j<h;++j){
            int y0 = y[i][j] - 1; // internal 1-based -> output 0-based
            fout << "X_"<<i<<"_"<<j<<" 0\n";
            fout << "Y_"<<i<<"_"<<j<<" "<< y0 <<"\n";
        }
    }
    fout.close();
    cout << "[Policy] wrote placement to: " << fn << "\n";
}

int main(int argc, char** argv){
    if(argc<4){
        cout<<"Usage: ./build/oc_policy_singlecol m h rounds "
              "[--T=12] [--policy=momentum|lagrange|irl1|soft] "
              "[--freeze=1] "
              "[--alpha=0.9 --eta=1.0 --top=0.10] "
              "[--rho=0.3] [--lam=2.0 --delta=2.0] [--gamma=0.3 --beta=3.0] "
              "[--cap=64 --scale=1000] [--dV=1] [--pow2=1] [--progress=1] "
              "[--homotopy=1 --sigma=0.8 --shrink=0.9] "
              "[--fast_ctr=1]\n";
        return 0;
    }
    // 物理 strip 尺寸（来自命令行）
    int m_full=atoi(argv[1]), h_full=atoi(argv[2]), rounds=atoi(argv[3]);
    long long T=12; string policy="soft";
    bool freeze=true; double alpha=0.9,eta=1.0,top=0.10,cap=64,scale=1000;
    double rho=0.3, lam=2.0, delta=2.0, gamma=0.3, beta=3.0;
    int dV=1; bool pow2=true, progress=true;
    bool homotopy=true; double sigma=0.8, shrink=0.9;
    bool fast_ctr=true; // m>h 时自动 shrink 到 h×h 再 Row-swap 扩展

    string v;
    for(int i=4;i<argc;++i){
        string s(argv[i]);
        if(parse_flag(s,"--T",v)) T=stoll(v);
        else if(parse_flag(s,"--policy",v)) policy=v;
        else if(parse_flag(s,"--freeze",v)) freeze=(v!="0");
        else if(parse_flag(s,"--alpha",v)) alpha=stod(v);
        else if(parse_flag(s,"--eta",v))   eta=stod(v);
        else if(parse_flag(s,"--top",v))   top=stod(v);
        else if(parse_flag(s,"--cap",v))   cap=stod(v);
        else if(parse_flag(s,"--scale",v)) scale=stod(v);
        else if(parse_flag(s,"--rho",v))   rho=stod(v);
        else if(parse_flag(s,"--lam",v))   lam=stod(v);
        else if(parse_flag(s,"--delta",v)) delta=stod(v);
        else if(parse_flag(s,"--gamma",v)) gamma=stod(v);
        else if(parse_flag(s,"--beta",v))  beta=stod(v);
        else if(parse_flag(s,"--dV",v))    dV=stoi(v);
        else if(parse_flag(s,"--pow2",v))  pow2=(v!="0");
        else if(parse_flag(s,"--progress",v)) progress=(v!="0");
        else if(parse_flag(s,"--homotopy",v)) homotopy=(v!="0");
        else if(parse_flag(s,"--sigma",v)) sigma=stod(v);
        else if(parse_flag(s,"--shrink",v)) shrink=stod(v);
        else if(parse_flag(s,"--fast_ctr",v)) fast_ctr=(v!="0");
    }

    // DP 用的尺寸（CTR 模式下缩成 h_full×h_full）
    int m_dp = m_full;
    int h_dp = h_full;
    bool use_ctr = fast_ctr && (m_full > h_full);
    if(use_ctr){
        m_dp = h_full;
        h_dp = h_full;
    }

    // Build policy（基于 DP 尺寸）
    WeighterParams P; P.export_scale=scale; P.cap_max=cap;
    P.alpha=alpha; P.eta=eta; P.top_ratio=top; P.rho=rho; P.lam=lam; P.delta=delta; P.gamma=gamma; P.beta=beta;

    IWeighterPolicy* core=nullptr;
    if(policy=="momentum") core=MakeMomentumPolicy(m_dp,h_dp,P);
    else if(policy=="lagrange") core=MakeLagrangePolicy(m_dp,h_dp,P);
    else if(policy=="irl1") core=MakeIRL1Policy(m_dp,h_dp,P);
    else if(policy=="soft") core=MakeSoftPolicy(m_dp,h_dp,P);
    else { cerr<<"Unknown policy: "<<policy<<"\n"; return 2; }
    if(freeze) core = MakeFreezeDecorator(core);

    // Initial weights = 1*scale  (DP 盘)
    Mat<long long> wH(m_dp, std::vector<long long>(h_dp-1, (long long)std::llround(scale)));
    Mat<long long> wV(m_dp-1, std::vector<long long>(h_dp,   (long long)std::llround(scale)));
    core->set_initial(wH,wV,(long long)std::llround(scale));

    cout<<"[Policy] m_full="<<m_full<<" h_full="<<h_full
        <<"  m_dp="<<m_dp<<" h_dp="<<h_dp
        <<"  rounds="<<rounds<<" T_target="<<T
        <<" policy="<<core->name()<<" cap="<<cap<<" scale="<<scale
        <<" freeze="<<(freeze?1:0)<<" homotopy="<<(homotopy?1:0)
        <<" sigma="<<sigma<<" shrink="<<shrink
        <<" fast_ctr="<<(fast_ctr?1:0)<<"\n";

    long long Tprime = -1; // 从 round1 的 maxL 冷启动
    Mat<int> y_best; long long best_hp=(long long)9e18;

    for(int it=1; it<=rounds; ++it){
        // 1) 导出权重并在 DP 盘上求解
        Mat<long long> WH,WV; core->export_ll(WH,WV);
        lightOCShortest::Config cfg;
        cfg.dV=(dV==0?1:dV); cfg.verbose=progress; cfg.progress=progress; cfg.max_steps=-1;
        cfg.omp_threads=0; cfg.low_mem=true; cfg.use_pow2_pack=pow2;
        cfg.W.enabled=true; cfg.W.wH=WH; cfg.W.wV=WV;

        lightOCShortest solver(cfg);
        auto t0=std::chrono::high_resolution_clock::now();
        auto R = solver.solve(m_dp,h_dp);
        auto t1=std::chrono::high_resolution_clock::now();
        double ms=std::chrono::duration<double,std::milli>(t1-t0).count();

        long long hpW = lightOCShortest::hpwl_sum_weighted(R.y_order, cfg.W, 1);
        long long hpE = lightOCShortest::hpwl_sum_equal   (R.y_order, 1);

        // 2) 在 DP 盘上计算 T 口径 maxL/WNS，用于 T′ 初始化/终止判断
        auto st_T = compute_maxL_WNS(R.y_order, m_dp, h_dp, T);
        if (homotopy && Tprime < 0){
            long long T0 = (long long)std::llround(std::ceil(sigma * (double)st_T.maxL));
            Tprime = std::max(T, T0);
            cout<<"[Init-Homotopy] round1 maxL="<<st_T.maxL
                <<"  set T_prime0="<<Tprime<<" (T_target="<<T<<")\n";
        }
        if (!homotopy && Tprime < 0) Tprime = T;

        // 3) 用当前 T′ 更新权重
        auto rep_Tprime = core->update_from_layout(R.y_order, Tprime, /*pick_by_ratio=*/true);

        // 4) 打印 T′ 和 T 两套口径
        auto st_Tprime = compute_maxL_WNS(R.y_order, m_dp, h_dp, Tprime);
        cout<<"  [Round "<<it<<"] T_prime="<<Tprime<<"  T_target="<<T
            <<"  HPWL_w="<<hpW<<"  HPWL_eq="<<hpE
            <<"  maxL="<<st_T.maxL
            <<"  WNS(T')="<<st_Tprime.WNS
            <<"  WNS*(T)="<<st_T.WNS
            <<"  picked="<<rep_Tprime.picked_edges
            <<"  solve_ms="<<ms<<"\n";

        if (hpW < best_hp){ best_hp=hpW; y_best=R.y_order; }

        // 5) 若对 T 已可行（maxL<=T），直接收敛
        if (st_T.WNS >= 0){
            core->freeze_if_feasible(true);
            y_best = R.y_order; best_hp = hpW;
            cout<<"  [TD] reach target on DP grid: all edges <= T_target.\n";
            break;
        }

        // 6) 仅当当前解已满足 T′ 时（maxL<=T′），才允许收紧 T′，并重置惯性
        if (homotopy && st_Tprime.WNS >= 0){
            long long maxL_curr  = st_T.maxL;
            long long cand_sigma  = (long long)std::llround(std::ceil(sigma  * (double)maxL_curr));
            long long cand_shrink = (long long)std::llround(std::ceil(shrink * (double)Tprime));
            long long T_next = std::max(T, std::min(cand_sigma, cand_shrink));

            if (T_next < Tprime){
                cout<<"    [Tighten] T_prime: "<<Tprime<<" -> "<<T_next
                    <<" (sigma*maxL="<<cand_sigma<<", shrink="<<cand_shrink<<")"
                    <<"  | cold-start: reset inertia\n";
                Tprime = T_next;
                core->reset_inertia();  // 冷启动：保留当前权重，清空动量
            }
        }
    }

    // 7) CTR 扩展：若启用 fast_ctr，则从 DP 盘（h×h）扩展到物理 strip（m_full×h_full）
    if(!y_best.empty()){
        Mat<int> y_out;
        if(use_ctr){
            y_out = RowSwapCTRExpander::expand(y_best, m_full, h_full);
        }else{
            y_out = y_best;
        }
        write_singlecol_placement(y_out, m_full, h_full, T);
    }
    delete core;
    return 0;
}
