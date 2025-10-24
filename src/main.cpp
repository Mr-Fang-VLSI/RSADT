#include "lightOCShortest.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <algorithm>
using namespace std;

// ===== RSAD seed (for UB print only; no pruning) =====
struct RSADSeed {
    static double L_I(int m,int h,int g){
        return -(2.0/3.0)*g*g*g + 2.0*h*g*g + ((2.0/3.0)-h*h-h)*g + m*1.0*h*h + m*h - m - h;
    }
    static int O_LL(int i,int j){ return (i>=j)? i*i - i + j : (j-1)*(j-1) + i; }
    static int O_LR(int i,int j,int g){
        if(i<=g-j+1) return int(-0.5*j*j + (g+1.5)*j - g + i - 1);
        return int(0.5*g*(g-1) + 0.5*i*(i-1) + j);
    }
    static int O_UL(int i,int j,int g){
        if(i<=g-j+1) return int(-0.5*i*i + (g+1.5)*i - g + j - 1);
        return int(0.5*g*(g-1) + 0.5*j*(j-1) + i);
    }
    static int O_UR(int i,int j,int g){
        if(i<=j) return 2*g*i - i*i + i + j - 2*g;
        else     return 2*g*j - j*j + i - g;
    }
    static int O_LM(int i,int j,int g){ return (j-1)*g + i; }
    static int O_CTR(int i,int j,int h){ return (i-1)*h + j; }

    struct Rec{int I,J,i,j,ord;};
    static vector<vector<int>> build(int m,int n,int g){
        vector<vector<int>> y(m, vector<int>(n,0));
        vector<Rec> LL,LM,LR,CTR,UL,UM,UR;
        auto push=[&](vector<Rec>& V,int I,int J,int i,int j,int O){ V.push_back({I,J,i,j,O}); };
        for(int I=0;I<g;++I){
            for(int J=0;J<g;++J)           push(LL,I,J, I+1,       J+1,           O_LL(I+1,J+1));
            for(int J=g;J<n-g;++J)         push(LM,I,J, I+1,      (J-g)+1,        O_LM(I+1,(J-g)+1,g));
            for(int J=n-g;J<n;++J)         push(LR,I,J, I+1,      (J-(n-g))+1,    O_LR(I+1,(J-(n-g))+1,g));
        }
        for(int I=g; I<m-g; ++I)
            for(int J=0; J<n; ++J)         push(CTR,I,J,(I-g)+1,  J+1,            O_CTR((I-g)+1,J+1,n));
        for(int I=m-g; I<m; ++I){
            for(int J=0;J<g;++J)           push(UL,I,J,(I-(m-g))+1, J+1,          O_UL((I-(m-g))+1,J+1,g));
            for(int J=g;J<n-g;++J)         push(UM,I,J,(I-(m-g))+1,(J-g)+1,       O_LM((I-(m-g))+1,(J-g)+1,g));
            for(int J=n-g;J<n;++J)         push(UR,I,J,(I-(m-g))+1,(J-(n-g))+1,   O_UR((I-(m-g))+1,(J-(n-g))+1,g));
        }
        auto sort_by=[&](vector<Rec>& V, const string& name){
            if(name=="CTR"){
                sort(V.begin(),V.end(),[](const Rec&a,const Rec&b){
                    if(a.ord!=b.ord) return a.ord<b.ord; if(a.i!=b.i) return a.i<b.i; return a.j<b.j; });
            }else if(name=="LM"||name=="UM"){
                sort(V.begin(),V.end(),[](const Rec&a,const Rec&b){
                    if(a.ord!=b.ord) return a.ord<b.ord; if(a.j!=b.j) return a.j<b.j; return a.i<b.i; });
            }else if(name=="LL"||name=="UL"){
                sort(V.begin(),V.end(),[](const Rec&a,const Rec&b){
                    if(a.ord!=b.ord) return a.ord<b.ord; if(a.i!=b.i) return a.i<b.i; return a.j<b.j; });
            }else{
                sort(V.begin(),V.end(),[](const Rec&a,const Rec&b){
                    if(a.ord!=b.ord) return a.ord<b.ord; if(a.j!=b.j) return a.j<b.j; return a.i<b.i; });
            }
        };
        sort_by(LL,"LL"); sort_by(LM,"LM"); sort_by(LR,"LR");
        sort_by(CTR,"CTR"); sort_by(UL,"UL"); sort_by(UM,"UM"); sort_by(UR,"UR");
        int cur=0;
        auto stamp=[&](vector<Rec>& V){ for(size_t k=0;k<V.size();++k){ y[V[k].I][V[k].J]=cur+(int)k+1; } cur+=(int)V.size(); };
        stamp(LL); stamp(LM); stamp(LR); stamp(CTR); stamp(UL); stamp(UM); stamp(UR);
        return y;
    }
    static vector<vector<int>> run(int m,int n,int& g_star, double& Ltheory){
        int gmax = std::max(1, std::min(m,n)/2);
        g_star = 1; Ltheory = L_I(m,n,1);
        for(int G=2; G<=gmax; ++G){ double v = L_I(m,n,G); if(v < Ltheory){ Ltheory=v; g_star=G; } }
        return build(m,n,g_star);
    }
};

static bool load_matrix(const string& path, int R, int C, vector<vector<long long>>& M){
    if(path.empty()) return false;
    ifstream fin(path);
    if(!fin) return false;
    M.assign(R, vector<long long>(C,0));
    long long x; int r=0, c=0;
    while(fin>>x){ M[r][c]=x; if(++c==C){ c=0; ++r; if(r==R) break; } }
    return (r==R && c==0);
}

static void run_case(int m,int h,long long dV,int steps,int threads,
                     bool use_pow2,bool check_oc,bool verify,
                     bool use_rsad,double ub_after_frac,
                     const string& wH_path,const string& wV_path)
{
    lightOCShortest::Config cfg;
    cfg.dV = dV; cfg.verbose = true; cfg.progress = true;
    cfg.max_steps = steps; cfg.omp_threads = threads;
    cfg.low_mem = true; cfg.use_pow2_pack = use_pow2;

    if(!wH_path.empty() && !wV_path.empty()){
        cfg.W.enabled = true;
        if(!load_matrix(wH_path, m, h-1, cfg.W.wH)){ cerr<<"[Warn] wH load failed; fallback equal\n"; cfg.W.enabled=false; }
        else if(!load_matrix(wV_path, m-1, h, cfg.W.wV)){ cerr<<"[Warn] wV load failed; fallback equal\n"; cfg.W.enabled=false; }
    }

    cout << "\n=== OC-Shortest m="<<m<<" h="<<h
         << " steps="<<(steps>0?steps:m*h)
         << " thr="<<(threads>0?threads:0)
         << " pow2="<<(use_pow2?1:0)
         << " weighted="<<(cfg.W.enabled?1:0)
         << " ===\n";

    if(use_rsad){
        int g_star=1; double Lth=0.0;
        auto y0 = RSADSeed::run(m,h,g_star,Lth);
        long long UB = cfg.W.enabled ? lightOCShortest::hpwl_sum_weighted(y0, cfg.W, 1)
                                     : lightOCShortest::hpwl_sum_equal(y0, 1);
        cout << "[RSAD] g*="<<g_star<<", L_theory="<< (long long)std::llround(Lth)
             << ", UB="<< UB << " (print only; no pruning)"
             << ", ub_after="<<ub_after_frac<<"\n";
    }

    lightOCShortest solver(cfg);
    auto t0 = std::chrono::high_resolution_clock::now();
    auto R  = solver.solve(m,h);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double,std::milli>(t1-t0).count();
    cout << "Solve-only time: " << ms << " ms\n";

    if(verify){
        long long HP = cfg.W.enabled ? lightOCShortest::hpwl_sum_weighted(R.y_order, cfg.W, 1)
                                     : lightOCShortest::hpwl_sum_equal(R.y_order, 1);
        cout << "[Verify] HPWL="<<HP<<" (solver_report="<<R.total_cost<<")\n";
    }
    if(check_oc){
        bool ok = lightOCShortest::check_OC_lin(R.y_order);
        cout << "[OC] " << (ok?"OK":"FAIL") << "\n";
    }
}

int main(int argc,char** argv){
    bool use_pow2=true, check_oc=true, verify=true, use_rsad=true;
    double ub_after=0.40;
    string wH_path, wV_path;

    if(argc<3){
        cout << "Usage: ./build/oc_shortest m h [dV=1] [steps=-1] [threads=0]\n"
                "  --pow2=1/0        use 2^k bit-pack key (default=1)\n"
                "  --check-oc=1/0    OC verification in main (default=1)\n"
                "  --verify=1/0      HPWL verify in main (default=1)\n"
                "  --rsad=1/0        print RSAD seed UB (default=1)\n"
                "  --ub-after=0.40   print-only suggestion (no pruning)\n"
                "  --wH path         H weights file: m x (h-1)\n"
                "  --wV path         V weights file: (m-1) x h\n";
        run_case(4,4,1,-1,0,use_pow2,check_oc,verify,use_rsad,ub_after,wH_path,wV_path);
        return 0;
    }

    int m = std::atoi(argv[1]);
    int h = std::atoi(argv[2]);
    long long dV = (argc>=4 ? std::atoll(argv[3]) : 1);
    int steps    = (argc>=5 ? std::atoi(argv[4]) : -1);
    int thr      = (argc>=6 ? std::atoi(argv[5]) : 0);

    for(int i=1;i<argc;++i){
        std::string s(argv[i]);
        auto eat=[&](const std::string& k,const std::string& v){
            if(k=="--pow2") use_pow2=(v=="1");
            else if(k=="--check-oc") check_oc=(v=="1");
            else if(k=="--verify") verify=(v=="1");
            else if(k=="--rsad") use_rsad=(v=="1");
            else if(k=="--ub-after") ub_after=stod(v);
        };
        auto pos=s.find('=');
        if(pos!=std::string::npos) eat(s.substr(0,pos), s.substr(pos+1));
        else if(s=="--wH" && i+1<argc) wH_path=argv[++i];
        else if(s=="--wV" && i+1<argc) wV_path=argv[++i];
    }

    run_case(m,h,dV,steps,thr,use_pow2,check_oc,verify,use_rsad,ub_after,wH_path,wV_path);
    return 0;
}
