#pragma once
#include <vector>

// RowSwapCTRExpander
// 用于把在“小盘”（通常是 h×h，如 8×8）上求出的 y 布局
// 扩展成在 tall strip（m_full×h_full）上的布局。
// 要求：只沿着行方向扩展，列数保持不变。
class RowSwapCTRExpander {
public:
    using MatI = std::vector<std::vector<int>>;

    // y_small: m_dp × h_dp 的 DP 布局（1-based label）
    // m_full, h_full: 目标条带尺寸（如 32×8, 40×8）
    //
    // 若形状不匹配（如 m_full<=m_dp 或 h_full!=h_dp），则直接返回 y_small。
    static MatI expand(const MatI& y_small, int m_full, int h_full);
};
