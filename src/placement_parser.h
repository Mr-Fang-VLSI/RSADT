#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>

class PlacementParser {
public:
    // —— 已有：单列导出（保留） ——
    static bool dump_single_column(
        const std::vector<std::vector<int>>& y_order, // m x h (1-based)
        int m, int h,
        const std::string& out_path,
        const std::string& header = std::string(),
        int x_index = 0,
        bool y_from_one = true)
    {
        if ((int)y_order.size()!=m || (m>0 && (int)y_order[0].size()!=h)) {
            std::cerr << "[Parser] dimension mismatch in single_column.\n";
            return false;
        }
        std::ofstream fout(out_path);
        if(!fout){ std::cerr<<"[Parser] cannot open: "<<out_path<<"\n"; return false; }
        if(!header.empty()) fout<<header<<"\n\n";
        for(int i=0;i<m;++i){
            for(int j=0;j<h;++j){
                int yval = y_order[i][j];
                if(y_from_one) yval -= 1;
                fout << "X_"<<i<<"_"<<j<<" " << x_index << "\n";
                fout << "Y_"<<i<<"_"<<j<<" " << yval    << "\n";
            }
        }
        return true;
    }

    // —— 新增：多列导出（X/Y 直接给到） ——
    static bool dump_multi_columns(
        const std::vector<std::vector<int>>& X, // m x h，列索引（0..K-1）
        const std::vector<std::vector<int>>& Y, // m x h，0-based 行号
        int m, int h,
        const std::string& out_path,
        const std::string& header = std::string())
    {
        if((int)X.size()!=m || (m>0 && (int)X[0].size()!=h)) { std::cerr<<"[Parser] X dim mismatch.\n"; return false; }
        if((int)Y.size()!=m || (m>0 && (int)Y[0].size()!=h)) { std::cerr<<"[Parser] Y dim mismatch.\n"; return false; }
        std::ofstream fout(out_path);
        if(!fout){ std::cerr<<"[Parser] cannot open: "<<out_path<<"\n"; return false; }
        if(!header.empty()) fout<<header<<"\n\n";
        for(int i=0;i<m;++i){
            for(int j=0;j<h;++j){
                fout << "X_"<<i<<"_"<<j<<" " << X[i][j] << "\n";
                fout << "Y_"<<i<<"_"<<j<<" " << Y[i][j] << "\n";
            }
        }
        return true;
    }

    // 文件名
    static std::string default_filename(int m,int h,long long T){
        std::ostringstream oss;
        if(m==h) oss<<"placement"<<m<<"_Col_1_T_"<<T<<".txt";
        else     oss<<"placement"<<m<<"x"<<h<<"_Col_1_T_"<<T<<".txt";
        return oss.str();
    }
    static std::string default_filename_cols(int m,int h,int cols,long long T){
        std::ostringstream oss;
        if(m==h) oss<<"placement"<<m<<"_Col_"<<cols<<"_T_"<<T<<".txt";
        else     oss<<"placement"<<m<<"x"<<h<<"_Col_"<<cols<<"_T_"<<T<<".txt";
        return oss.str();
    }

    // 头部描述
    static std::string default_header(int m,int h,int cols=1){
        std::ostringstream oss;
        oss << "This is a sample " << m << "x" << h
            << " PE array being placed on " << cols
            << " DSP " << (cols>1?"columns:":"column:");
        return oss.str();
    }
};
