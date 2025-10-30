#include "strip_placer.h"
#include "placement_parser.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <cctype>

using std::cout; using std::cerr; using std::endl; using std::string;
template<typename T> using Mat = std::vector<std::vector<T>>;

static bool parse_flag(const string& s, const string& key, string& out){
    auto p=s.find('='); if(p==string::npos) return false;
    if(s.substr(0,p)==key){ out=s.substr(p+1); return true; }
    return false;
}
static bool is_integer(const char* s){
    if(!s||!*s) return false; const char* p=s;
    if(*p=='+'||*p=='-') ++p; if(!*p) return false;
    while(*p){ if(!std::isdigit((unsigned char)*p)) return false; ++p; }
    return true;
}

int main(int argc,char** argv){
    if(argc<4){
        cout<<"Usage: ./build/oc_split m h K [rounds=5] [dV=1]\n"
              "  --T=20      Max-T for strip solver\n"
              "  --alpha=0.5 --eta=0.5 --top=0.10 --cap=32 --scale=1000\n"
              "  --pow2=1    use 2^k bitpack (default 1)\n"
              "  --progress=1\n"
              "  --outfile path\n";
        return 0;
    }
    const int m = std::atoi(argv[1]);
    const int h = std::atoi(argv[2]);
    int   K     = std::atoi(argv[3]);
    int rounds  = (argc>=5 && is_integer(argv[4]) ? std::atoi(argv[4]) : 5);
    long long dV= (argc>=6 && is_integer(argv[5]) ? std::atoll(argv[5]) : 1);

    StripPlacer::Params P; P.m=m; P.h=h; P.K=K; P.rounds=rounds; P.dV=dV;
    P.T = std::max(m, h/K);
    string out_path; string v;
    for(int i=1;i<argc;++i){
        string s(argv[i]);
        if(parse_flag(s,"--T",v))       P.T=std::llround(std::stod(v));
        else if(parse_flag(s,"--alpha",v)) P.alpha=std::stod(v);
        else if(parse_flag(s,"--eta",v))   P.eta  =std::stod(v);
        else if(parse_flag(s,"--top",v))   P.top  =std::stod(v);
        else if(parse_flag(s,"--cap",v))   P.cap  =std::stod(v);
        else if(parse_flag(s,"--scale",v)) P.scale=std::stod(v);
        else if(parse_flag(s,"--pow2",v))  P.use_pow2=(v!="0");
        else if(parse_flag(s,"--progress",v)) P.progress=(v=="1");
        else if(s=="--outfile" && i+1<argc) out_path=argv[++i];
    }

    cout<<"[Split] m="<<m<<" h="<<h<<" K="<<K<<" w="<<(h/K)
        <<" rounds="<<P.rounds<<" dV="<<P.dV<<" T="<<P.T
        <<" pow2="<<(P.use_pow2?1:0)<<endl;

    try{
        StripPlacer placer(P);
        placer.solve_once();
        Mat<int> X, Y; placer.build_global_xy(X, Y);

        if(out_path.empty()){
            out_path = PlacementParser::default_filename_cols(m,h,K,P.T);
        }
        string header = PlacementParser::default_header(m,h,K);
        if(!PlacementParser::dump_multi_columns(X, Y, m, h, out_path, header)){
            cerr<<"[Split] write output failed.\n"; return 2;
        }
        cout<<"[Split] wrote placement to: "<<out_path<<endl;
    }catch(const std::exception& e){
        cerr<<"[Split] ERROR: "<<e.what()<<endl;
        return 3;
    }
    return 0;
}
