#include "multicol_splitter.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <fstream>
#include <sstream>
#include <iostream>
#include <regex>
#include <algorithm>

using std::string;

bool MultiColSplitter::ensure_dir(const std::string& dir){
#ifdef _WIN32
    int rc = _mkdir(dir.c_str());
#else
    int rc = mkdir(dir.c_str(), 0755);
#endif
    if (rc == 0) return true;
    if (errno == EEXIST) return true;
    std::perror(("mkdir "+dir).c_str());
    return false;
}

int MultiColSplitter::run_cmd(const std::string& cmd){
    std::cout << "[MultiCol] exec: " << cmd << std::endl;
    int rc = std::system(cmd.c_str());
    if (rc == -1){
        std::perror("system");
        return rc;
    }
#ifdef _WIN32
    return rc;
#else
    if (WIFEXITED(rc)) return WEXITSTATUS(rc);
    return rc;
#endif
}

bool MultiColSplitter::parse_singlecol_placement(const std::string& path,
                                                 int m, int w,
                                                 std::vector<std::vector<long long>>& y_base)
{
    std::ifstream fin(path);
    if(!fin){
        std::cerr << "[MultiCol] ERROR: cannot open base placement: " << path << std::endl;
        return false;
    }
    y_base.assign(m, std::vector<long long>(w, -1));

    // Lines like: "Y_i_j <y>"
    std::regex reY(R"(Y_(\d+)_(\d+)\s+(-?\d+))");
    std::string line;
    long long cnt = 0;
    while(std::getline(fin, line)){
        std::smatch mth;
        if(std::regex_match(line, mth, reY) && mth.size()==4){
            int i = std::stoi(mth[1].str());
            int j = std::stoi(mth[2].str());
            long long y = std::stoll(mth[3].str());
            if(i>=0 && i<m && j>=0 && j<w){
                y_base[i][j] = y;
                cnt++;
            }
        }
    }
    // Check filled
    for(int i=0;i<m;++i){
        for(int j=0;j<w;++j){
            if(y_base[i][j] < 0){
                std::cerr << "[MultiCol] ERROR: missing Y_"<<i<<"_"<<j<<" in "<<path<<std::endl;
                return false;
            }
        }
    }
    std::cout << "[MultiCol] parsed base Y entries: " << cnt << " (expect "<<(m*w)<<")" << std::endl;
    return true;
}

bool MultiColSplitter::write_multicol_placement(const std::string& out_path,
                                                int m, int h, int K, int w,
                                                const std::vector<std::vector<long long>>& y_base,
                                                int snake)
{
    std::ofstream fout(out_path);
    if(!fout){
        std::cerr << "[MultiCol] ERROR: cannot write " << out_path << std::endl;
        return false;
    }
    fout << "This is a sample " << m << "x" << h
         << " PE array being placed on " << K << " DSP columns"
         << (snake? " (snake=ON: odd strips use logical L-R flip)":" (snake=OFF)") << ":\n\n";

    for(int s=0; s<K; ++s){
        const bool flip_lr = (snake != 0) && (s % 2 == 1);
        for(int i=0;i<m;++i){
            for(int j_local=0;j_local<w;++j_local){
                const int j_global = s*w + j_local;
                const int j_idx = flip_lr ? (w - 1 - j_local) : j_local; // logical left-right flip
                const long long y = y_base[i][j_idx];

                fout << "X_"<<i<<"_"<<j_global<<" " << s << "\n";
                fout << "Y_"<<i<<"_"<<j_global<<" " << y << "\n";
            }
        }
    }
    std::cout << "[MultiCol] wrote multi-col placement: " << out_path << std::endl;
    return true;
}

std::string MultiColSplitter::run(const MultiColOptions& opt)
{
    // 这里原来只允许 K=2/4，这里放宽到 2/4/5，其它 K 仍视为非法
    if (opt.m<=0 || opt.h<=0 || (opt.K!=2 && opt.K!=4 && opt.K!=5)){
        std::cerr << "[MultiCol] ERROR: invalid m/h/K. m="<<opt.m<<" h="<<opt.h<<" K="<<opt.K<<std::endl;
        return "";
    }
    if (opt.h % opt.K != 0){
        std::cerr << "[MultiCol] ERROR: h="<<opt.h<<" cannot be evenly split by K="<<opt.K<<std::endl;
        return "";
    }
    const int w = opt.h / opt.K;
    const int m = opt.m;
    const int T = opt.T;

    // working dir for base run
    std::string work = opt.work_dir.empty() ? "." : opt.work_dir;
    if (!ensure_dir(work)) return "";

    // --- Compose command for single-col run on m x w ---
    std::ostringstream oss;
    oss << opt.single_bin
        << " " << m << " " << w << " " << opt.rounds
        << " --T=" << T
        << " --policy=" << opt.policy
        << " --freeze=" << opt.freeze
        << " --cap=" << opt.cap
        << " --scale=" << opt.scale
        << " --dV=" << opt.dV
        << " --pow2=" << opt.pow2
        << " --progress=" << opt.progress;

    if (opt.policy == "soft"){
        oss << " --alpha="<<opt.alpha
            << " --eta="<<opt.eta
            << " --gamma="<<opt.gamma
            << " --beta="<<opt.beta;
        if (opt.homotopy) {
            oss << " --homotopy=1"
                << " --sigma="<<opt.sigma
                << " --shrink="<<opt.shrink;
        }
    } else if (opt.policy == "lagrange"){
        oss << " --rho="<<opt.rho;
        if (opt.homotopy){
            oss << " --homotopy=1"
                << " --sigma="<<opt.sigma
                << " --shrink="<<opt.shrink;
        }
    } else if (opt.policy == "irl1"){
        oss << " --lam="<<opt.lam
            << " --delta="<<opt.delta;
        if (opt.homotopy){
            oss << " --homotopy=1"
                << " --sigma="<<opt.sigma
                << " --shrink="<<opt.shrink;
        }
    } else {
        std::cerr << "[MultiCol] WARN: unknown policy="<<opt.policy<<", fallback run anyway.\n";
    }

    int rc = run_cmd(oss.str());
    if (rc != 0){
        std::cerr << "[MultiCol] ERROR: single-col solver failed, rc="<<rc<<std::endl;
        return "";
    }

    // move base placement into work_dir (for parsing)
    std::ostringstream base_name;
    base_name << "placement" << m << "_Col_1_T_" << T << ".txt";
    std::string base_src = base_name.str();

    std::ostringstream base_dst_oss;
    base_dst_oss << work << "/base_placement_"<< m << "x" << w << "_T_"<<T<<".txt";
    std::string base_dst = base_dst_oss.str();

    {
        std::ifstream in(base_src, std::ios::binary);
        if(!in){
            std::cerr << "[MultiCol] ERROR: cannot find base placement: " << base_src << std::endl;
            return "";
        }
        std::ofstream out(base_dst, std::ios::binary);
        out << in.rdbuf();
    }
    std::remove(base_src.c_str());
    std::cout << "[MultiCol] moved base placement to: " << base_dst << std::endl;

    // parse base Y map
    std::vector<std::vector<long long>> y_base;
    if (!parse_singlecol_placement(base_dst, m, w, y_base)){
        return "";
    }

    // write replicated multi-col placement (with optional logical L-R flip on odd strips)
    std::ostringstream final_out;
    final_out << "placement" << m << "_Col_"<< opt.K << "_T_" << T << ".txt";
    std::string out_path = final_out.str();
    if (!write_multicol_placement(out_path, m, opt.h, opt.K, w, y_base, opt.snake)){
        return "";
    }
    return out_path;
}
