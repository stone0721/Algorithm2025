#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <functional> // 提供 std::hash 用于哈希函数
#include <fstream>
using namespace std;

// 定义 pair_hash 结构，用于 std::unordered_map 的键（std::pair）哈希
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        // 对 pair 的两个元素分别计算哈希值，并通过异或运算组合
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// 定义 Segment 结构，表示查询序列和参考序列之间的匹配片段信息
struct Segment {
    int32_t x0;       // 片段在查询序列中的起始 x 坐标
    int32_t y0;       // 片段在参考序列中的起始 y 坐标
    int32_t length;   // 片段的长度
    int8_t direction; // 方向：1 表示主对角线（直接匹配），0 表示反对角线（互补匹配）
    int32_t distance; // 片段中的错误数（不匹配的位点数）
};

// 全局变量，存储点阵图和对角线匹配表
std::vector<std::vector<uint8_t>> grid; // 点阵图，记录查询序列和参考序列的匹配状态
std::vector<std::vector<int16_t>> main_ld, anti_lu; // 主对角线和反对角线累计匹配表

// 检查两个碱基是否互补（A-T 或 C-G）
bool is_complement(char a, char b) {
    return (a == 'A' && b == 'T') || (a == 'T' && b == 'A') ||
           (a == 'C' && b == 'G') || (a == 'G' && b == 'C');
}

// 构建点阵图，比较查询序列和参考序列
// 参数：query（查询序列），reference（参考序列）
// 返回：M x N 的点阵图，1 表示直接匹配，2 表示互补匹配，0 表示不匹配
std::vector<std::vector<uint8_t>> build_dotplot(const std::string& query, const std::string& reference) {
    int M = query.length(), N = reference.length(); // 获取序列长度
    // 初始化 M x N 的点阵图，初始值为 0
    std::vector<std::vector<uint8_t>> result(M, std::vector<uint8_t>(N, 0));
    
    for (int x = 0; x < M; ++x) {
        for (int y = 0; y < N; ++y) {
            if (query[x] == reference[y]) {
                result[x][y] = 1; // 碱基相同，标记为直接匹配
            } else if (is_complement(query[x], reference[y])) {
                result[x][y] = 2; // 碱基互补，标记为互补匹配
            }
        }
    }
    return result;
}

// 初始化主对角线和反对角线匹配表
// 参数：grid（点阵图）
// 返回：包含主对角线表（main_ld）和反对角线表（anti_lu）的 pair
std::pair<std::vector<std::vector<int16_t>>, std::vector<std::vector<int16_t>>>
init_diag_tables(const std::vector<std::vector<uint8_t>>& grid) {
    int M = grid.size(), N = grid[0].size(); // 获取点阵图尺寸
    std::vector<std::vector<int16_t>> main_ld(M, std::vector<int16_t>(N, 0)); // 主对角线表
    std::vector<std::vector<int16_t>> anti_lu(M, std::vector<int16_t>(N, 0)); // 反对角线表

    // 主对角线：从左下到右上，记录连续直接匹配的累计长度
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            int16_t cnt = (grid[i][j] == 1) ? 1 : 0; // 当前点是否为直接匹配
            if (i > 0 && j > 0) cnt += main_ld[i-1][j-1]; // 加上前一个对角线点的累计值
            main_ld[i][j] = cnt;
        }
    }

    // 反对角线：从右上到左下，记录连续互补匹配的累计长度
    for (int i = 0; i < M; ++i) {
        for (int j = N-1; j >= 0; --j) {
            int16_t cnt = (grid[i][j] == 2) ? 1 : 0; // 当前点是否为互补匹配
            if (i > 0 && j + 1 < N) cnt += anti_lu[i-1][j+1]; // 加上前一个反对角线点的累计值
            anti_lu[i][j] = cnt;
        }
    }
    return {main_ld, anti_lu};
}

// 寻找对角线上的匹配片段
// 参数：grid（点阵图）
// 返回：包含所有匹配片段的向量
std::vector<Segment> find_diagonal_segments_np(const std::vector<std::vector<uint8_t>>& grid) {
    int M = grid.size(), N = grid[0].size(); // 获取点阵图尺寸
    std::vector<Segment> segs; // 存储匹配片段

    // 辅助函数：从对角线提取连续匹配段
    auto extract_runs = [](const std::vector<uint8_t>& line) {
        std::vector<int> starts, ends; // 记录匹配段的起点和终点
        std::vector<int8_t> mask(line.size() + 2, 0); // 创建掩码，标记匹配点
        for (size_t i = 0; i < line.size(); ++i) mask[i + 1] = (line[i] != 0); // 非零表示匹配
        for (size_t i = 1; i < mask.size() - 1; ++i) {
            if (mask[i] == 1 && mask[i-1] == 0) starts.push_back(i-1); // 匹配段起点
            if (mask[i] == 1 && mask[i+1] == 0) ends.push_back(i); // 匹配段终点
        }
        return std::make_pair(starts, ends);
    };

    // 主对角线
    for (int d = -M + 1; d < N; ++d) { // 遍历所有可能的主对角线，d = y - x
        std::vector<uint8_t> diag; // 存储对角线上的网格值
        int start_x = std::max(0, -d), start_y = start_x + d; // 计算起始坐标
        for (int x = start_x, y = start_y; x < M && y < N; ++x, ++y) {
            diag.push_back(grid[x][y]); // 收集对角线上的值
        }
        auto [starts, ends] = extract_runs(diag); // 提取连续匹配段
        for (size_t i = 0; i < starts.size(); ++i) {
            int length = ends[i] - starts[i]; // 计算匹配段长度
            if (length > 0) { // 仅保留有效长度的片段
                int x0 = start_x + starts[i]; // 片段起点 x 坐标
                int y0 = x0 + d; // 片段起点 y 坐标
                segs.push_back({x0, y0, length, 1, 0}); // 添加片段
            }
        }
    }

    // 反对角线（通过水平翻转网格处理）
    std::vector<std::vector<uint8_t>> flipped(M, std::vector<uint8_t>(N)); // 翻转网格
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            flipped[i][j] = grid[i][N-1-j]; // 水平翻转
        }
    }
    for (int d = -M + 1; d < N; ++d) {
        std::vector<uint8_t> diag; // 存储反对角线上的网格值
        int start_x = std::max(0, -d), start_y = N - 1 - (start_x + d); // 计算起始坐标
        for (int x = start_x, y = start_y; x < M && y >= 0; ++x, --y) {
            diag.push_back(flipped[x][N-1-y]); // 收集翻转后的网格值
        }
        auto [starts, ends] = extract_runs(diag); // 提取连续匹配段
        for (size_t i = 0; i < starts.size(); ++i) {
            int length = ends[i] - starts[i]; // 计算匹配段长度
            if (length > 0) {
                int x0 = start_x + starts[i]; // 片段起点 x 坐标
                int y0 = N - 1 - (start_x + starts[i] + d); // 片段起点 y 坐标
                segs.push_back({x0, y0, length, 0, 0}); // 添加片段
            }
        }
    }
    return segs;
}

// 加权区间调度算法，选择不重叠的片段以最大化总长度
// 参数：segs（输入片段向量）
// 返回：最优片段集合
std::vector<Segment> weighted_interval_scheduling_np(const std::vector<Segment>& segs) {
    if (segs.empty()) return {}; // 空输入返回空向量

    // 按片段结束位置（x1）排序
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 + segs[a].length - 1 < segs[b].x0 + segs[b].length - 1;
    });

    // 存储片段的起点、终点和长度
    std::vector<int> x0s(segs.size()), x1s(segs.size()), lens(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        x0s[i] = segs[order[i]].x0;
        lens[i] = segs[order[i]].length;
        x1s[i] = x0s[i] + lens[i] - 1;
    }

    // 计算每个片段的前驱（最后一个不重叠的片段）
    std::vector<int> p(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        auto it = std::lower_bound(x1s.begin(), x1s.end(), x0s[i]);
        p[i] = (it == x1s.begin()) ? -1 : std::distance(x1s.begin(), it) - 1;
    }

    // 动态规划选择最优片段
    std::vector<int> dp(segs.size() + 1, 0); // dp[j] 表示前 j 个片段的最大总长度
    std::vector<bool> choose(segs.size(), false); // 标记是否选择该片段
    for (size_t j = 1; j <= segs.size(); ++j) {
        int w = lens[j-1]; // 当前片段的长度
        if (w + dp[p[j-1] + 1] > dp[j-1]) { // 选择当前片段是否更优
            dp[j] = w + dp[p[j-1] + 1];
            choose[j-1] = true;
        } else {
            dp[j] = dp[j-1];
        }
    }

    // 回溯构建结果
    std::vector<Segment> result;
    int j = segs.size();
    while (j > 0) {
        if (choose[j-1]) {
            result.push_back(segs[order[j-1]]); // 选择该片段
            j = p[j-1] + 1;
        } else {
            --j;
        }
    }
    std::reverse(result.begin(), result.end()); // 反转以按起始位置升序
    return result;
}

// 在片段之间的空白区域合并片段
// 参数：segs（输入片段），rate（最大错误率）
// 返回：合并后的片段集合
std::vector<Segment> merge_in_blanks_np(const std::vector<Segment>& segs, float rate) {
    if (segs.empty()) return {}; // 空输入返回空向量

    // 按 x0 升序排序
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 < segs[b].x0;
    });

    std::vector<Segment> merged; // 存储合并结果
    auto prev = segs[order[0]]; // 前一个片段
    for (size_t i = 1; i < segs.size(); ++i) {
        auto s = segs[order[i]]; // 当前片段
        if (s.direction == prev.direction) { // 确保方向相同
            // 检查是否在同一对角线
            bool same_diag = (s.direction == 1 && (s.y0 - s.x0) == (prev.y0 - prev.x0)) ||
                             (s.direction == 0 && (s.x0 + s.y0) == (prev.x0 + prev.y0));
            if (same_diag) {
                int new_len = (s.x0 + s.length - 1) - prev.x0 + 1; // 合并后的长度
                int prev_count, affective; // 计算匹配数
                if (prev.direction == 1) { // 主对角线
                    prev_count = (prev.x0 > 0 && prev.y0 > 0) ? main_ld[prev.x0-1][prev.y0-1] : 0;
                    int sx1 = s.x0 + s.length - 1;
                    int sy1 = s.y0 + (s.length - 1);
                    affective = main_ld[sx1][sy1] - prev_count; // 使用累计表计算匹配数
                } else { // 反对角线
                    prev_count = (prev.x0 > 0 && prev.y0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[prev.x0-1][prev.y0+1] : 0;
                    int sx1 = s.x0 + s.length - 1;
                    int sy1 = s.y0 - (s.length - 1);
                    affective = anti_lu[sx1][sy1] - prev_count;
                }
                int distance = new_len - affective; // 错误数
                if (distance / static_cast<float>(new_len) < rate) { // 错误率符合要求
                    prev.length = new_len;
                    prev.distance = distance;
                    continue; // 继续尝试合并下一个片段
                }
            }
        }
        merged.push_back(prev); // 保存前一个片段
        prev = s; // 更新前一个片段
    }
    merged.push_back(prev); // 添加最后一个片段
    return merged;
}

// 在指定范围内寻找满足错误率要求的大片段
// 参数：x_start, x_end（范围），rate（最大错误率），min_len（最小长度）
// 返回：找到的片段集合
std::vector<Segment> find_large_valid_segments_in_range_np(int x_start, int x_end, float rate, int min_len) {
    int M = grid.size(), N = grid[0].size(); // 获取点阵图尺寸
    std::vector<Segment> recs; // 存储结果
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int>>, pair_hash> groups; // 按对角线分组

    // 主对角线
    for (int d = -(M-1); d < N; ++d) { // 遍历所有可能的主对角线
        int lo = std::max({0, -d, x_start}); // 范围下界
        int hi = std::min({M-1, N-1-d, x_end}); // 范围上界
        if (hi - lo + 1 >= 1) { // 确保范围有效
            std::vector<std::pair<int, int>> pts; // 存储对角线上的点
            for (int x = lo; x <= hi; ++x) {
                pts.emplace_back(x, x + d); // 添加 (x, y) 坐标
            }
            groups[std::make_pair(1, d)] = pts; // 按方向和对角线编号存储
        }
    }

    // 反对角线
    for (int d = 0; d < M + N - 1; ++d) { // 遍历所有可能的反对角线
        int lo = std::max({0, d - (N-1), x_start});
        int hi = std::min({M-1, d, x_end});
        if (hi - lo + 1 >= 1) {
            std::vector<std::pair<int, int>> pts;
            for (int x = lo; x <= hi; ++x) {
                pts.emplace_back(x, d - x);
            }
            groups[std::make_pair(0, d)] = pts;
        }
    }

    // 处理每条对角线
    for (const auto& group : groups) {
        int direction = group.first.first; // 对角线方向
        const auto& pts = group.second; // 对角线上的点
        int L = pts.size(); // 点数
        std::vector<int> mis(L + 1, 0); // 累计错误数
        for (int i = 0; i < L; ++i) {
            auto [x, y] = pts[i];
            mis[i+1] = mis[i] + ((direction == 1 && grid[x][y] != 1) || (direction == 0 && grid[x][y] != 2)); // 记录错误
        }

        // 动态规划寻找最优片段
        std::vector<int> dp(L + 1, 0); // dp[R] 表示前 R 个点的最大长度
        std::vector<int> prev_idx(L + 1, 0); // 记录前驱索引
        int L_ptr = 0; // 滑动窗口左端
        for (int R = 1; R <= L; ++R) {
            while (L_ptr < R && (mis[R] - mis[L_ptr]) > rate * (R - L_ptr)) { // 确保错误率
                ++L_ptr;
            }
            dp[R] = dp[R-1];
            prev_idx[R] = R-1;
            if (R - L_ptr >= min_len) { // 长度满足要求
                int length = R - L_ptr;
                if (dp[L_ptr] + length > dp[R]) { // 更新最优解
                    dp[R] = dp[L_ptr] + length;
                    prev_idx[R] = L_ptr;
                }
            }
        }

        // 回溯提取片段
        int idx = L;
        while (idx > 0) {
            int pi = prev_idx[idx];
            if (pi < idx - 1) {
                int start = pi, end = idx - 1;
                auto [x0, y0] = pts[start]; // 片段起点
                int length = end - start + 1; // 片段长度
                int errs = mis[end + 1] - mis[start]; // 错误数
                recs.push_back(Segment{x0, y0, length, static_cast<int8_t>(direction), errs});
                idx = pi;
            } else {
                --idx;
            }
        }
    }
    return recs;
}

// 全局填充空白区域
// 参数：segs（输入片段），rate（最大错误率），min_gap（最小间隙长度）
// 返回：填充后的片段集合
std::vector<Segment> fill_in_blanks_global_np(const std::vector<Segment>& segs, float rate, int min_gap) {
    if (segs.empty()) return {};

    // 按 x0 升序排序
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 < segs[b].x0;
    });

    std::vector<Segment> recs; // 存储结果
    for (size_t i = 0; i < segs.size() - 1; ++i) {
        auto prev = segs[order[i]]; // 前一个片段
        auto curr = segs[order[i + 1]]; // 当前片段
        recs.push_back(prev);
        int g0 = prev.x0 + prev.length; // 前一片段结束
        int g1 = curr.x0 - 1; // 当前片段开始前
        int gap_len = g1 - g0 + 1; // 间隙长度
        if (gap_len >= min_gap) { // 间隙足够大
            auto extras = find_large_valid_segments_in_range_np(g0, g1, rate, min_gap); // 寻找间隙中的片段
            recs.insert(recs.end(), extras.begin(), extras.end()); // 添加新片段
        }
    }
    recs.push_back(segs[order.back()]); // 添加最后一个片段

    recs = merge_in_blanks_np(recs, rate); // 合并片段
    recs = weighted_interval_scheduling_np(recs); // 优化选择
    return recs;
}

// 向后扩展片段的结束位置
// 参数：segs（输入片段），rate（最大错误率）
// 返回：扩展后的片段集合
std::vector<Segment> extend_end_backward_np(const std::vector<Segment>& segs, float rate) {
    if (segs.empty()) return {};

    // 按 x0 降序排序
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 > segs[b].x0;
    });

    std::vector<Segment> out; // 存储结果
    int prev_end = segs[order[0]].x0 + segs[order[0]].length; // 前一片段结束位置
    for (size_t i = 0; i < segs.size(); ++i) {
        auto seg = segs[order[i]];
        int ox0 = seg.x0, oy0 = seg.y0, olen = seg.length, odir = seg.direction, odist = seg.distance;
        int ox1 = ox0 + olen - 1; // 片段结束 x 坐标
        int oy1 = oy0 + (olen - 1) * (odir == 1 ? 1 : -1); // 片段结束 y 坐标

        int target_end = prev_end - 1; // 目标结束位置
        int space = target_end - ox1; // 扩展空间
        int new_len = olen, new_dist = odist; // 初始化新长度和错误数
        int cand_x1 = ox1, cand_y1 = oy1; // 候选结束坐标

        if (ox1 < target_end) { // 可以向后扩展
            cand_x1 = target_end;
            int step = (odir == 1) ? 1 : -1; // 根据方向调整 y 坐标
            cand_y1 = oy1 + step * space;
            new_len = olen + space;

            while (cand_x1 > ox1) { // 尝试扩展
                if (cand_x1 >= static_cast<int>(grid.size()) || cand_y1 < 0 || cand_y1 >= static_cast<int>(grid[0].size())) { // 边界检查
                    --cand_x1;
                    cand_y1 -= step;
                    --new_len;
                    continue;
                }
                int prev_corr, corr; // 计算匹配数
                if (odir == 1) {
                    prev_corr = (ox0 > 0 && oy0 > 0) ? main_ld[ox0-1][oy0-1] : 0;
                    corr = main_ld[cand_x1][cand_y1] - prev_corr;
                } else {
                    prev_corr = (ox0 > 0 && oy0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[ox0-1][oy0+1] : 0;
                    corr = anti_lu[cand_x1][cand_y1] - prev_corr;
                }
                int dist = new_len - corr; // 错误数
                if (dist / static_cast<float>(new_len) < rate) { // 错误率符合要求
                    new_dist = dist;
                    break;
                }
                --cand_x1;
                cand_y1 -= step;
                --new_len;
            }
        }
        out.push_back(Segment{ox0, oy0, new_len, static_cast<int8_t>(odir), new_dist});
        prev_end = ox0; // 更新前一片段结束位置
    }

    // 按 x0 升序重新排序
    std::vector<size_t> asc_order(out.size());
    for (size_t i = 0; i < out.size(); ++i) asc_order[i] = i;
    std::sort(asc_order.begin(), asc_order.end(), [&out](size_t a, size_t b) {
        return out[a].x0 < out[b].x0;
    });
    std::vector<Segment> result(out.size());
    for (size_t i = 0; i < out.size(); ++i) result[i] = out[asc_order[i]];
    return result;
}

// 向前扩展片段的起始位置
// 参数：segs（输入片段），rate（最大错误率）
// 返回：扩展后的片段集合
std::vector<Segment> extend_start_backward_np(const std::vector<Segment>& segs, float rate) {
    if (segs.empty()) return {};

    // 按 x0 升序排序
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 < segs[b].x0;
    });

    std::vector<Segment> out; // 存储结果
    int prev_end = -1; // 前一片段结束位置
    for (size_t i = 0; i < segs.size(); ++i) {
        auto seg = segs[order[i]];
        int ox0 = seg.x0, oy0 = seg.y0, olen = seg.length, odir = seg.direction, odist = seg.distance;
        int ox1 = ox0 + olen - 1;
        int oy1 = oy0 + (olen - 1) * (odir == 1 ? 1 : -1);

        int target_start = prev_end + 1; // 目标起始位置
        int space = ox0 - target_start; // 扩展空间
        int cand_x0 = ox0, cand_y0 = oy0, cand_len = olen; // 候选起始坐标和长度

        if (ox0 > target_start) { // 可以向前扩展
            cand_x0 = target_start;
            int step = (odir == 1) ? 1 : -1;
            cand_y0 = oy0 - step * space;
            cand_len = olen + space;

            while (cand_x0 < ox0) { // 尝试扩展
                if (cand_x0 < 0 || cand_y0 < 0 || cand_y0 >= static_cast<int>(grid[0].size())) break; // 边界检查
                int prev_corr, corr;
                if (odir == 1) {
                    prev_corr = (cand_x0 > 0 && cand_y0 > 0) ? main_ld[cand_x0-1][cand_y0-1] : 0;
                    corr = main_ld[cand_x0 + cand_len - 1][cand_y0 + cand_len - 1] - prev_corr;
                } else {
                    prev_corr = (cand_x0 > 0 && cand_y0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[cand_x0-1][cand_y0+1] : 0;
                    corr = anti_lu[cand_x0 + cand_len - 1][cand_y0 - (cand_len - 1)] - prev_corr;
                }
                int dist = cand_len - corr;
                if (dist / static_cast<float>(cand_len) < rate) { // 错误率符合要求
                    out.push_back(Segment{cand_x0, cand_y0, cand_len, static_cast<int8_t>(odir), dist});
                    break;
                }
                ++cand_x0;
                cand_y0 += step;
                --cand_len;
            }
            if (cand_x0 == ox0) out.push_back(seg); // 未成功扩展，使用原片段
        } else {
            out.push_back(seg); // 无需扩展，直接添加
        }
        prev_end = out.back().x0 + out.back().length - 1; // 更新前一片段结束位置
    }
    return out;
}

// 选择长度大于等于 min_len 的片段
// 参数：segs（输入片段），min_len（最小长度）
// 返回：满足条件的片段集合
std::vector<Segment> chose_segs_np(const std::vector<Segment>& segs, int min_len) {
    std::vector<Segment> result;
    for (const auto& seg : segs) {
        if (seg.length >= min_len) result.push_back(seg); // 仅保留长度足够的片段
    }
    return result;
}

// 最小区间覆盖算法，选择覆盖查询序列的最少片段
// 参数：segs（输入片段），rate（最大错误率），length_thresh（最小长度）
// 返回：覆盖片段集合
std::vector<Segment> minimal_interval_cover2_np(const std::vector<Segment>& segs, float rate, int length_thresh) {
    if (segs.empty()) return {};

    // 确定覆盖范围
    int min_x = segs[0].x0, max_x = segs[0].x0 + segs[0].length - 1;
    for (const auto& s : segs) {
        min_x = std::min(min_x, s.x0);
        max_x = std::max(max_x, s.x0 + s.length - 1);
    }

    // 构建桶，记录每个起点的最大结束位置
    std::vector<int> best_end(max_x + 1, -1);
    std::vector<int> best_idx(max_x + 1, -1);
    for (size_t idx = 0; idx < segs.size(); ++idx) {
        int bx0 = segs[idx].x0;
        int bx1 = bx0 + segs[idx].length - 1;
        if (bx1 > best_end[bx0]) { // 更新最大结束位置
            best_end[bx0] = bx1;
            best_idx[bx0] = idx;
        }
    }

    std::vector<Segment> result; // 存储结果
    int covered_end = min_x - 1; // 已覆盖的结束位置
    while (covered_end < max_x) {
        int start_pos = covered_end + 1; // 下一次覆盖的起点
        int candidate_end = -1, candidate_i = -1; // 候选结束位置和索引
        for (int x0 = min_x; x0 <= start_pos; ++x0) {
            if (best_end[x0] > candidate_end) { // 寻找最远结束位置
                candidate_end = best_end[x0];
                candidate_i = best_idx[x0];
            }
        }

        if (candidate_i < 0 || candidate_end < start_pos) break; // 无法继续覆盖

        auto best = segs[candidate_i]; // 最佳片段
        int new_len = std::min(best.x0 + best.length - 1, candidate_end) - start_pos + 1; // 新片段长度
        int new_y0, prev, correct; // 计算匹配数
        if (best.direction == 1) {
            new_y0 = best.y0 + (start_pos - best.x0); // 新 y 坐标
            prev = (start_pos > 0 && new_y0 > 0) ? main_ld[start_pos-1][new_y0-1] : 0;
            correct = main_ld[start_pos + new_len - 1][new_y0 + new_len - 1] - prev;
        } else {
            new_y0 = best.y0 - (start_pos - best.x0);
            prev = (start_pos > 0 && new_y0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[start_pos-1][new_y0+1] : 0;
            correct = anti_lu[start_pos + new_len - 1][new_y0 - new_len + 1] - prev;
        }
        int distance = new_len - correct; // 错误数

        if (new_len >= length_thresh && distance / static_cast<float>(new_len) < rate) { // 满足长度和错误率要求
            result.push_back(Segment{start_pos, new_y0, new_len, static_cast<int8_t>(best.direction), distance});
        }
        covered_end = candidate_end; // 更新覆盖结束位置
    }
    return result;
}

// 合并具有容差的片段
// 参数：segs（输入片段），max_gap（最大间隙），max_error_rate（最大错误率）
// 返回：合并后的片段集合
std::vector<Segment> merge_with_tolerance_np(const std::vector<Segment>& segs, int max_gap, float max_error_rate) {
    if (segs.empty()) return {};

    // 按方向和对角线编号分组
    std::unordered_map<std::pair<int, int>, std::vector<size_t>, pair_hash> groups;
    std::vector<int> diag_id(segs.size());
    std::vector<int> x1(segs.size()), y1(segs.size());

    for (size_t i = 0; i < segs.size(); ++i) {
        x1[i] = segs[i].x0 + segs[i].length - 1; // 片段结束 x 坐标
        y1[i] = segs[i].direction == 1 ? segs[i].y0 + segs[i].length - 1 : segs[i].y0 - (segs[i].length - 1); // 结束 y 坐标
        diag_id[i] = segs[i].direction == 1 ? segs[i].y0 - segs[i].x0 : segs[i].y0 + segs[i].x0; // 对角线编号
        groups[std::make_pair(segs[i].direction, diag_id[i])].push_back(i); // 分组
    }

    std::vector<Segment> merged; // 存储合并结果
    for (const auto& group : groups) {
        int dir_val = group.first.first; // 对角线方向
        auto indices = group.second; // 片段索引

        // 按 x0 排序
        std::sort(indices.begin(), indices.end(), [&segs](size_t a, size_t b) {
            return segs[a].x0 < segs[b].x0;
        });

        auto curr = segs[indices[0]]; // 当前片段
        int cx1 = curr.x0 + curr.length - 1; // 当前片段结束 x 坐标
        int cy1 = curr.direction == 1 ? curr.y0 + curr.length - 1 : curr.y0 - (curr.length - 1); // 结束 y 坐标

        for (size_t i = 1; i < indices.size(); ++i) {
            auto s = segs[indices[i]];
            int sx1 = s.x0 + s.length - 1; // 新片段结束 x 坐标
            int sy1 = s.direction == 1 ? s.y0 + s.length - 1 : s.y0 - (s.length - 1); // 结束 y 坐标
            int gap = s.x0 - cx1 - 1; // 间隙长度
            int merged_len = sx1 - curr.x0 + 1; // 合并后长度
            bool cond_align = (dir_val == 1 && s.y0 - cy1 == gap + 1) ||
                              (dir_val == 0 && cy1 - s.y0 == gap + 1); // 检查对角线对齐

            if (gap <= max_gap && (gap + curr.distance + s.distance) / static_cast<float>(merged_len) <= max_error_rate && cond_align) {
                // 合并片段
                curr.length = merged_len;
                curr.distance = gap + curr.distance + s.distance;
                cx1 = curr.x0 + curr.length - 1;
                cy1 = curr.y0 + (curr.length - 1) * (curr.direction == 1 ? 1 : -1);
            } else {
                merged.push_back(curr); // 保存当前片段
                curr = s; // 更新当前片段
                cx1 = curr.x0 + curr.length - 1;
                cy1 = curr.y0 + (curr.length - 1) * (curr.direction == 1 ? 1 : -1);
            }
        }
        merged.push_back(curr); // 添加最后一个片段
    }
    return merged;
}

// 打印匹配结果
// 参数：matches（匹配片段集合）
// 输出：长度 >= 30 的片段坐标信息
void get_answer_np(const std::vector<Segment>& matches) {
    std::cout << "Remain segments:" << std::endl;
    for (const auto& seg : matches) {
        if (seg.length >= 30) { // 仅打印长度 >= 30 的片段
            int x0 = seg.x0, length = seg.length, y0 = seg.y0;
            int x1 = x0 + length - 1; // 结束 x 坐标
            if (seg.direction == 1) { // 主对角线
                int y1 = y0 + length - 1; // 结束 y 坐标
                std::cout << "(" << x0 << "," << x1 + 1 << "," << y0 << "," << y1 + 1 << "),";
            } else { // 反对角线
                int y1 = y0 - (length - 1);
                std::cout << "(" << x0 << "," << x1 + 1 << "," << y1 << "," << y0 + 1 << "),";
            }
        }
    }
    std::cout << std::endl;
}

// 打印详细输出，计算最终得分
// 参数：matches（匹配片段集合）
// 输出：每个片段的详细信息和总得分
void get_detail_np(const std::vector<Segment>& matches) {
    int point = 0;
    for (const auto& seg : matches) {
        if (seg.length <= 0) continue;
        int x0 = seg.x0, y0 = seg.y0, length = seg.length, direction = seg.direction;
        int x1 = x0 + length - 1;
        int affective;
        if (direction == 1) {
            int y1 = y0 + length - 1;
            int prev = (x0 > 0 && y0 > 0) ? main_ld[x0-1][y0-1] : 0;
            affective = main_ld[x1][y1] - prev;
            std::cout << "(" << x0 << "," << x1 + 1 << "," << y0 << "," << y1 + 1 << ")\t";
        } else {
            int y1 = y0 - (length - 1);
            int prev = (x0 > 0 && y0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[x0-1][y0+1] : 0;
            affective = anti_lu[x1][y1] - prev;
            std::cout << "(" << x0 << "," << x1 + 1 << "," << y1 << "," << y0 + 1 << ")\t";
        }
        std::cout << ";length = " << length << "\taffective length = " << affective << "\t";
        if (affective / static_cast<float>(length) > 0.9) {
            std::cout << "True" << std::endl;
            point += affective;
        } else {
            std::cout << "False" << std::endl;
        }
    }
    std::cout << "final score = " << point << std::endl;
}



// 主函数，寻找最佳匹配片段
// 参数：query（查询序列），reference（参考序列）
// 返回：最佳匹配片段集合
std::vector<Segment> find_best_matches_np(const std::string& query, const std::string& reference) {
    grid = build_dotplot(query, reference); // 构建点阵图
    auto [ml, al] = init_diag_tables(grid); // 初始化主对角线和反对角线匹配表
    main_ld = ml;
    anti_lu = al;

    // 执行一系列处理步骤
    auto final_segs = find_diagonal_segments_np(grid); // 提取对角线上的匹配段
    final_segs = merge_with_tolerance_np(final_segs, 1, 0.039); // 合并具有小间隙的段
    final_segs = minimal_interval_cover2_np(final_segs, 0.08, 20); // 最小区间覆盖
    final_segs = fill_in_blanks_global_np(final_segs, 0.065, 25); // 全局填充空白区域
    final_segs = chose_segs_np(final_segs, 25); // 选择长度 >= 25 的段
    final_segs = extend_start_backward_np(final_segs, 0.1); // 向前扩展起始

    // 按 x0 降序排序以进行向后扩展
    std::vector<size_t> order(final_segs.size());
    for (size_t i = 0; i < final_segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&final_segs](size_t a, size_t b) {
        return final_segs[a].x0 > final_segs[b].x0;
    });
    std::vector<Segment> temp(final_segs.size());
    for (size_t i = 0; i < final_segs.size(); ++i) temp[i] = final_segs[order[i]];
    final_segs = extend_end_backward_np(temp, 0.1); // 向后扩展结束

    // 按 x0 升序重新排序
    std::sort(final_segs.begin(), final_segs.end(), [](const Segment& a, const Segment& b) {
        return a.x0 < b.x0;
    });

    return final_segs;
}

int main() {
    std::string reference = "TGATTTAGAACGGACTTAGCAGACATTGAAACTCGAGGGGTATAGCAATAGATGCCCAAAAAGGTAAGCGCCATAAGCGTGGTTCTACGAGCCAGGTGCTCATGCCTAAGTTCTGCGCCTTCGCTGTCACTTGGAAATACTGTAATGGATCATGCCTAGGTTATGCGCCTTCGGGGTCACTTCAACATACTGTAATGGATCATGCCTAGGTTTTGCGTGTTCGCTGTCATTTCGAAATACTCCAATGGATGATGCCTAGGTTCTGTGCCTTCGCTGACGCATGGAAATACTTTAACGGATCATGCCCAGGCTCTGCGCCTTCGCTGAAACTTCGAAATACTCTAATGGATCATGCCTCGGTGCTCCACCTTCGCTTTCATTCCGAAATACTCTAATGGATCGCGTCCGTGTAACAACTTCGTACTGTTATATCGACCATCAGAATACCCATCCCTCGGGGAGGTAACCTATATTCACGTCGCAAGTTTCGATCTACAGTATGCTGACTGTTTGCCGCGATTTTAAGTCAAGAAGCGCGTCCACATGGTCGCGGTCGTCAACTTCAGTACGCTCATATGACACCAAAAGATCTACCTACAGCCCGTGCAGCTCGACTTTTGTGCTCTAGGGCACGACGGGTGGCGTTTGCTCCCGCGCATCTCGACTTTTAAGCTCTATGGCACAACGTGTGGCGTTTGCCCCCGCGCAGCTCGACTTTTGTGCTCTAGGGCACGGCGGGTGGCGTTTGCCCTCGCCCAGCTTGACTTTTGTGCTCTAGGGCACGACGGGTGGCGTTTGCCCCCGTGCAGCCCGACTTTTGTACTCTAGTGCACGACGGGTGGCGTTTGCCCCCGCACCGCTCGACTTTTGTGATCTAGGGCACTACGAGTAGCGTTGGCCCAGACAGATCAACGCACATGGATTCTCTAACAGGCCCCGCGCTTCTCATTGGCCCGTGAGACGGGTCTGAGAGGAAGACATTAGGTAGATACGGAAAGCTTTGGCGTAGTTCGTATCTTTCAGCGGTGAAGCGTCTTCGGTCCGGGCTGCGTTATGCCTGCGGGAGGAAGGCTCCACTAGATGGTTTACGAGACATAATGTCAGCTTCCTAAAGGTACGGCAGCGCCTGCGTATATCACAGGACGAATTGTCAGTTTGCTAGGGGTACGGGAGCGCTTGCGTATTACATAGGACGAATCGTCAGCTTCCTAAAGGGACGGTAGCGCTTGCGTGTTACATAGGACGAATTGTCAGCTTCGTAAAGGTACGGTAGTTCTTGCGTATTACATAGGATGCATTGTCCGCTTCCTAAAGGTACGCTGGCGCTTGCGTATCACATAGGACGGATAGCGCGATTGCTAAAGGTACGGGAGCGCTTGCGTCTTAGAGCGCACGAATCGGATATAAGCTGCGCCCGCGTCTGGCGAGCAAAAATCGTGGGAGCCAGCGAGGGAAAAACTGCTCGGGCGACTTAAACGGAATTACAAGACTCATTGCCATCGAGGACGTTAGACTAAAGAGCCCCTGCGTGCCTCCTTTGTATAGCTCGATGTAGTGGCCCGTGTATGTGGAACAGGAATGCTCGATCTAAGGTAGTAGTGGCTACAGCTCCGAGAGTTTGCGTACTGCGGTGCCAGGGATTTTGCCTGCGGGTACAGCCTCTGCGCACGCCGGTCTGTGATCAAGAACTAAACTAGAGA";
    std::string query = "TGATTTAGAACGGACTTAGCAGACATTGAAACTCGAGGGGTATAGCAATAGATGCCCAAAAAGGTAAGCGCCATAAGCGTGTTTCTACGAGCCAGGTGCTCATGCCTAAGTTCTGCGCCTTCGCTGTCACTGGGAAATACTGTAATGGATCATCCGTAGGTTATGCGCCTTCGGGGTCACTTCAACATACTGTAACGGATCGTGCCTAGGTTTTGCGTATTCGCTGTCATTTCGAATTACACCAATGGATGATGCCTAGGTTCTGTGCCTCCGCTGACGCATCGAAATACTTTAACGGATCGCGTCCGAGTAACAACTTCGTACTGTTATATAGGCAATCAGAATACCCATGCCTCGGGGAGGTAACCTATATTCACGTCGCAAGTTTCGATCTACAGTACTGTAGGTATATCTTTTGGTGTCATATGAGGGTACTGAACTTGACGACCGCGACCATGTGGATGCGCTTCTTGACTTAAAATCGCGGCAAACAGTAAGCATCCGTGAAGCTCGACTTTTGTGCTCTAGGGCACGACGGGTGGCGTTTGCTCCCGCGCATCTCGAGTTGTAAGCTCTATGGCACAACGGGTGGCGTTTGCCGCCGAGCAGCTCGACTTTTGTGCTCTAGGGCACGGCGGGTGGCGTTTGCCCTCGCCCAGCTTGACTTTTGTGCTCTAGGGCACGACGGGTGGCCTTTGCCCCCGCGCAGCTGGACTTTTGTGCTCTAGGGCACGGCGGGTGGCGTTTGCCCTCCCCCAGCTTGACTACTGTGCTCTAGGGCACGACGGGTGGCGTTTGCCCCCGCGCAGCTCGACTTTTGTGCTCTATGGCACGGGGGGTGGCGTTTGCCCTCGCCCAGCTTGACTTTTGCGCTCTAGGGCACGACGGGTGGCGTTTGCCGGCAAACGCCACACGTCGTGCCCTAGAGCACAAACGTCAAGCTGGGCGAGGGCAACCGCCACCCGCCCTGCCCTAGAGCACAAAAGTCGAGCTGCGCGGGCCCGCGCAGCTCGACTTTTGTGCTCTAGGACACGGCGGGTGGCGTTTGCCCTCGCCCAGCTTGACTCTTGTGCTCTAGGGCACGACGGGTGGCGTTTGCCCCAGCGCAGCCCGACTTTTGTACTCTAGAGCACGACGGTTGGCATTTGCCCCCGCACCGCTCGACTTTTGTGATCTAGGGCCCTAGGAGTAGCGTTGGCCAGCTTTCCGTATCTACCTAATGTCTTCCTCTCAGACCCGTCTCACGGGCCAATGAGAAGCGCGGGGCCTATTAGAGAATCCATGTGCGTTGATCTGTCTGCAGACAGCTCAACGCACATGGATTCGCTAGCAGGCCCCGCGCTTCTCATTGGCCCGTGAGACGGGTCTGAGAGGAAGACATAAGGTAGATACGGCAAGCTCACGTCCGTGTAACAACGTCGTACTGTTATATCGACCATCAGAATCCCCATCCCGCGAGGAGGTAACCTATATTCAGGTCGCAAGTTTCGATCTACAGTATTGGCGTAGTTCGTATCTTTCAGCGGTGAAGCTTCTTCGGTCCGGGCTGCGTTATGCCTGCTGGAGGACGGCTCCACTAGATGGTTTACGAGACATAATGTCCGCTTCCTAAAGGTACACTGGCGCTTGAGTATCACATAGGACGGATAGCTCGATTCCTAAAGGGACGGGAGCGCTTGCGTCTTAGAGCGCATGAATCGTCAGCTTCCCAAAGGGACCGTAGCGCTTGCGTGTTATATAGGAAGAATGGTCAGCTTTGTAAAGGTACGGTAGTTCTTGCGTATTACAGAGGATGCATTGTCTACTACCTAAAGGTACGGCAGCGCCTGCGTATATCACAGGACGAATTGTCAGTTTGCTAGGGGTACGGGAGCGCTTGCATATTACATAGGACGAATCGGATATAAGCTGCGCCCGCGTCTGGCGATAAAAAATCGTGGTAGCCAGCGAGGGAAAAACTGCTCGGGCGACTTAAACGGAATTAAAAGACTCATTGCCGTGACAGACTTCCGTATAGCAACCTCTGGGATGTCGATGCGGTGTCCCCAGTCTGCGCTGAGCGGGGGCAGACAGACTTAGTTATAGTATGCATCTGTTTTAGCTAGACATCACGACCTAGTGGGGTTCATGTTGAGATTCTAGGGCGGTACGCAGCCGGTGGATTATTACTTCCCCAGAAATTCTGACTTCGTCACTGGATGGATTGTACTATCCGGTCAACCTTACAAGGTTTCAACAGGGACGAAGGGTAAACGTATGAAGCTTGGATGCCGTTACCGTAAAGGGCCCTATTGAAGTGTCGAGGACGTTAGACTAAAGAGCCCCTGCGTGCCTCCTTTGTATAGCTCGAGGTAGTGGCCCGGATATGTGGAACAGGAATGCTCGATCTAAGGTAGTAGTGGGTACCGCTCCGAGAGTTTGCGTACTGCGGTGCCCGGGATTTTGCCTGCGGGTACAGCCTCTGCGCACGCCGGTCTGTAATCAAGAACTAAACTAGAGA";
    auto matches = find_best_matches_np(query, reference);
    get_answer_np(matches);
    get_detail_np(matches);

    return 0;
}
