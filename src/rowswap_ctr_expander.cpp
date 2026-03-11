#include "rowswap_ctr_expander.h"
#include <algorithm>
#include <limits>

using MatI = std::vector<std::vector<int>>;

namespace {

// 8×8 基盘方向规范化：
// 期望 canonical 方向：上 4 行包含 label 0..31（即内部 1..32）。
// 若上 4 行最大值 !=31，则简单做一次转置修正。
void canonicalize_8x8(MatI& y) {
    const int m = (int)y.size();
    if (m == 0) return;
    const int h = (int)y[0].size();
    if (m != 8 || h != 8) return;

    int max_top0 = std::numeric_limits<int>::min();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 8; ++j) {
            int v0 = y[i][j] - 1; // 内部 1-based 转成 0-based
            if (v0 > max_top0) max_top0 = v0;
        }
    }

    // canonical 方向：上半区最大值为 31（0-based）
    if (max_top0 == 31) return;

    // 否则认为方向“歪了”，做一次转置
    MatI t(8, std::vector<int>(8, 0));
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            t[j][i] = y[i][j];
        }
    }
    y.swap(t);
}

struct CellRank {
    int i;
    int j;
    int val; // 1-based DP label
};

} // anonymous namespace

RowSwapCTRExpander::MatI
RowSwapCTRExpander::expand(const MatI& y_small_in, int m_full, int h_full)
{
    const int m_dp = (int)y_small_in.size();
    if (m_dp == 0) return y_small_in;
    const int h_dp = (int)y_small_in[0].size();

    // 只沿行方向扩展，列数必须一致；否则不做 CTR 拓展。
    if (h_dp != h_full || m_full <= m_dp) {
        return y_small_in;
    }

    MatI y_small = y_small_in; // 本地拷贝，用于方向修正等

    // 8×8 情况：先把方向标准化到“上半区是 0..31”
    canonicalize_8x8(y_small);

    const int m = m_full;
    const int h = h_full;

    // 从小盘的 top/bottom band 映射到大盘的 top/bottom band，中间是 CTR 行。
    const int top_src_rows    = m_dp / 2;          // 8 -> 4
    const int bottom_src_rows = m_dp - top_src_rows;
    const int top_big_rows    = top_src_rows;
    const int bottom_big_rows = bottom_src_rows;
    const int ctr_big_rows    = m - top_big_rows - bottom_big_rows;

    if (ctr_big_rows < 0) {
        // m_full 太小，无法插入 CTR 区，直接返回原解
        return y_small_in;
    }

    MatI y_big(m, std::vector<int>(h, 0));

    const long long n_full = (long long)m * h;
    const long long n_top  = (long long)top_big_rows    * h;
    const long long n_bot  = (long long)bottom_big_rows * h;
    const long long n_ctr  = n_full - n_top - n_bot;

    // -------- 收集 base square 的上下 band --------
    std::vector<CellRank> top_cells;
    std::vector<CellRank> bot_cells;
    top_cells.reserve((size_t)n_top);
    bot_cells.reserve((size_t)n_bot);

    for (int i = 0; i < top_src_rows; ++i) {
        for (int j = 0; j < h_dp; ++j) {
            top_cells.push_back({i, j, y_small[i][j]});
        }
    }
    const int start_bot_src = m_dp - bottom_src_rows;
    for (int i = start_bot_src; i < m_dp; ++i) {
        for (int j = 0; j < h_dp; ++j) {
            bot_cells.push_back({i, j, y_small[i][j]});
        }
    }

    auto cmp_val = [](const CellRank& a, const CellRank& b){ return a.val < b.val; };
    std::sort(top_cells.begin(), top_cells.end(), cmp_val);
    std::sort(bot_cells.begin(), bot_cells.end(), cmp_val);

    // -------- 顶部 band：1 .. n_top --------
    long long cur = 1;
    for (size_t k = 0; k < top_cells.size(); ++k, ++cur) {
        int src_i = top_cells[k].i; // 0..top_src_rows-1
        int src_j = top_cells[k].j; // 0..h-1
        int dst_i = src_i;          // 映射到最上面的 top_big_rows 行
        if (dst_i < 0 || dst_i >= m) continue;
        y_big[dst_i][src_j] = (int)cur;
    }

    // -------- 中间 CTR band：n_top+1 .. n_top+n_ctr --------
    for (int i = top_big_rows; i < top_big_rows + ctr_big_rows; ++i) {
        for (int j = 0; j < h; ++j) {
            if (cur <= n_top + n_ctr) {
                y_big[i][j] = (int)cur;
                ++cur;
            }
        }
    }

    // -------- 底部 band：最后 n_bot 个 label --------
    for (size_t k = 0; k < bot_cells.size() && cur <= n_full; ++k, ++cur) {
        int src_i = bot_cells[k].i;    // start_bot_src .. m_dp-1
        int src_j = bot_cells[k].j;

        int offset = src_i - start_bot_src;          // 0..bottom_src_rows-1
        int dst_i  = m - bottom_big_rows + offset;   // 映射到最底部 bottom_big_rows 行
        if (dst_i < 0 || dst_i >= m) continue;
        y_big[dst_i][src_j] = (int)cur;
    }

    return y_big;
}
