#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <functional>

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
    int16_t x0;       // 片段在查询序列中的起始 x 坐标
    int16_t y0;       // 片段在参考序列中的起始 y 坐标
    int16_t length;   // 片段的长度
    int8_t direction; // 方向：1 表示主对角线（直接匹配），0 表示反对角线（互补匹配）
    int16_t distance; // 片段中的错误数（不匹配的位点数）
};

// 声明全局变量，存储查询序列和参考序列的指针，以便在函数中高效访问
const std::string* query_ptr = nullptr;
const std::string* reference_ptr = nullptr;

// 动态计算网格中指定位置 (x, y) 的值
// 返回值：1 表示直接匹配（相同碱基），2 表示互补匹配（A-T 或 C-G），0 表示不匹配
uint8_t get_grid_value(int x, int y) {
    char q = (*query_ptr)[x], r = (*reference_ptr)[y]; // 获取查询序列和参考序列的字符
    if (q == r) return 1; // 直接匹配，例如 A-A、C-C
    if ((q == 'A' && r == 'T') || (q == 'T' && r == 'A') ||
        (q == 'C' && r == 'G') || (q == 'G' && r == 'C')) return 2; // 互补匹配
    return 0; // 无匹配
}

// 计算指定对角线范围内的匹配数
// 参数：x_start, y_start（起始坐标），length（对角线长度），direction（方向），M（查询序列长度），N（参考序列长度）
// 返回：匹配点的数量
int compute_diagonal_matches(int x_start, int y_start, int length, int direction, int M, int N) {
    int matches = 0; // 匹配计数
    for (int i = 0; i < length; ++i) {
        int x = x_start + i; // 查询序列坐标
        int y = direction == 1 ? y_start + i : y_start - i; // 参考序列坐标，主对角线 y 递增，反对角线 y 递减
        if (x >= M || y < 0 || y >= N) break; // 边界检查
        uint8_t val = get_grid_value(x, y); // 获取网格值
        matches += (direction == 1 && val == 1) || (direction == 0 && val == 2); // 累加匹配数
    }
    return matches;
}

// 寻找对角线上的匹配片段，不存储整个网格以节省内存
// 参数：M（查询序列长度），N（参考序列长度），min_len（最小片段长度，默认为5）
// 返回：包含所有匹配片段的向量
std::vector<Segment> find_diagonal_segments_np(int M, int N, int min_len = 5) {
    std::vector<Segment> segs; // 存储匹配片段

    // 处理主对角线（direction = 1，直接匹配）
    for (int d = -M + 1; d < N; ++d) { // 遍历所有可能的主对角线，d = y - x
        int start_x = std::max(0, -d), start_y = start_x + d; // 计算对角线的起始坐标
        std::vector<uint8_t> diag; // 存储对角线上的网格值
        for (int x = start_x, y = start_y; x < M && y < N; ++x, ++y) { // 沿对角线遍历
            diag.push_back(get_grid_value(x, y)); // 记录网格值
        }
        // 提取连续匹配段（run）
        std::vector<int> starts, ends; // 记录匹配段的起点和终点
        std::vector<int8_t> mask(diag.size() + 2, 0); // 创建掩码，标记匹配点（1 表示匹配）
        for (size_t i = 0; i < diag.size(); ++i) mask[i + 1] = (diag[i] == 1); // 填充掩码
        for (size_t i = 1; i < mask.size() - 1; ++i) {
            if (mask[i] == 1 && mask[i-1] == 0) starts.push_back(i-1); // 匹配段起点
            if (mask[i] == 1 && mask[i+1] == 0) ends.push_back(i); // 匹配段终点
        }
        for (size_t i = 0; i < starts.size(); ++i) {
            int length = ends[i] - starts[i]; // 计算匹配段长度
            if (length >= min_len) { // 仅保留长度达到阈值的片段
                int x0 = start_x + starts[i]; // 片段起点 x 坐标
                int y0 = x0 + d; // 片段起点 y 坐标
                segs.push_back({static_cast<int16_t>(x0), static_cast<int16_t>(y0), 
                               static_cast<int16_t>(length), 1, 0}); // 添加片段
            }
        }
    }

    // 处理反对角线（direction = 0，互补匹配）
    for (int d = -M + 1; d < N; ++d) { // 遍历所有可能的反对角线
        int start_x = std::max(0, -d), start_y = N - 1 - (start_x + d); // 计算起始坐标
        std::vector<uint8_t> diag; // 存储反对角线上的网格值
        for (int x = start_x, y = start_y; x < M && y >= 0; ++x, --y) { // 沿反对角线遍历
            diag.push_back(get_grid_value(x, N-1-y)); // 记录网格值（注意 y 坐标翻转）
        }
        // 提取连续匹配段
        std::vector<int> starts, ends;
        std::vector<int8_t> mask(diag.size() + 2, 0); // 创建掩码
        for (size_t i = 0; i < diag.size(); ++i) mask[i + 1] = (diag[i] == 2); // 标记互补匹配
        for (size_t i = 1; i < mask.size() - 1; ++i) {
            if (mask[i] == 1 && mask[i-1] == 0) starts.push_back(i-1); // 匹配段起点
            if (mask[i] == 1 && mask[i+1] == 0) ends.push_back(i); // 匹配段终点
        }
        for (size_t i = 0; i < starts.size(); ++i) {
            int length = ends[i] - starts[i]; // 计算匹配段长度
            if (length >= min_len) { // 仅保留长度达到阈值的片段
                int x0 = start_x + starts[i]; // 片段起点 x 坐标
                int y0 = N - 1 - (start_x + starts[i] + d); // 片段起点 y 坐标
                segs.push_back({static_cast<int16_t>(x0), static_cast<int16_t>(y0), 
                               static_cast<int16_t>(length), 0, 0}); // 添加片段
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

    // 按片段结束位置排序
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
// 参数：segs（输入片段），rate（最大错误率），M，N（序列长度）
// 返回：合并后的片段集合
std::vector<Segment> merge_in_blanks_np(const std::vector<Segment>& segs, float rate, int M, int N) {
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
                int affective = compute_diagonal_matches(prev.x0, prev.y0, new_len, prev.direction, M, N); // 匹配数
                int distance = new_len - affective; // 错误数
                if (distance / static_cast<float>(new_len) < rate) { // 错误率符合要求
                    prev.length = static_cast<int16_t>(new_len);
                    prev.distance = static_cast<int16_t>(distance);
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
// 参数：x_start, x_end（范围），rate（最大错误率），min_len（最小长度），M，N（序列长度）
// 返回：找到的片段集合
std::vector<Segment> find_large_valid_segments_in_range_np(int x_start, int x_end, float rate, int min_len, int M, int N) {
    std::vector<Segment> recs; // 存储结果
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int>>, pair_hash> groups; // 按对角线分组

    // 主对角线
    for (int d = -(M-1); d < N; ++d) { // 遍历所有可能的主对角线
        int lo = std::max({0, -d, x_start}); // 范围下界
        int hi = std::min({M-1, N-1-d, x_end}); // 范围上界
        if (hi - lo + 1 >= min_len) { // 确保范围足够长
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
        if (hi - lo + 1 >= min_len) {
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
            mis[i+1] = mis[i] + ((direction == 1 && get_grid_value(x, y) != 1) || (direction == 0 && get_grid_value(x, y) != 2)); // 记录错误
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
                recs.push_back({static_cast<int16_t>(x0), static_cast<int16_t>(y0), 
                               static_cast<int16_t>(length), static_cast<int8_t>(direction), 
                               static_cast<int16_t>(errs)});
                idx = pi;
            } else {
                --idx;
            }
        }
    }
    return recs;
}

// 全局填充空白区域
// 参数：segs（输入片段），rate（最大错误率），min_gap（最小间隙长度），M，N（序列长度）
// 返回：填充后的片段集合
std::vector<Segment> fill_in_blanks_global_np(const std::vector<Segment>& segs, float rate, int min_gap, int M, int N) {
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
            auto extras = find_large_valid_segments_in_range_np(g0, g1, rate, min_gap, M, N); // 寻找间隙中的片段
            recs.insert(recs.end(), extras.begin(), extras.end()); // 添加新片段
        }
    }
    recs.push_back(segs[order.back()]); // 添加最后一个片段

    recs = merge_in_blanks_np(recs, rate, M, N); // 合并片段
    recs = weighted_interval_scheduling_np(recs); // 优化选择
    return recs;
}

// 向后扩展片段的结束位置
// 参数：segs（输入片段），rate（最大错误率），M，N（序列长度）
// 返回：扩展后的片段集合
std::vector<Segment> extend_end_backward_np(const std::vector<Segment>& segs, float rate, int M, int N) {
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
                if (cand_x1 >= M || cand_y1 < 0 || cand_y1 >= N) { // 边界检查
                    --cand_x1;
                    cand_y1 -= step;
                    --new_len;
                    continue;
                }
                int corr = compute_diagonal_matches(ox0, oy0, new_len, odir, M, N); // 计算匹配数
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
        out.push_back({static_cast<int16_t>(ox0), static_cast<int16_t>(oy0), 
                      static_cast<int16_t>(new_len), static_cast<int8_t>(odir), 
                      static_cast<int16_t>(new_dist)});
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
// 参数：segs（输入片段），rate（最大错误率），M，N（序列长度）
// 返回：扩展后的片段集合
std::vector<Segment> extend_start_backward_np(const std::vector<Segment>& segs, float rate, int M, int N) {
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
                if (cand_x0 < 0 || cand_y0 < 0 || cand_y0 >= N) break; // 边界检查
                int corr = compute_diagonal_matches(cand_x0, cand_y0, cand_len, odir, M, N);
                int dist = cand_len - corr;
                if (dist / static_cast<float>(cand_len) < rate) { // 错误率符合要求
                    out.push_back({static_cast<int16_t>(cand_x0), static_cast<int16_t>(cand_y0), 
                                  static_cast<int16_t>(cand_len), static_cast<int8_t>(odir), 
                                  static_cast<int16_t>(dist)});
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
// 参数：segs（输入片段），rate（最大错误率），length_thresh（最小长度），M，N（序列长度）
// 返回：覆盖片段集合
std::vector<Segment> minimal_interval_cover2_np(const std::vector<Segment>& segs, float rate, int length_thresh, int M, int N) {
    if (segs.empty()) return {};

    // 确定覆盖范围
    int min_x = segs[0].x0, max_x = segs[0].x0 + segs[0].length - 1;
    for (const auto& s : segs) {
        min_x = std::min(min_x, static_cast<int>(s.x0));
        max_x = std::max(max_x, static_cast<int>(s.x0 + s.length - 1));
    }

    // 记录每个起点的最大结束位置
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
        int new_len = std::min(static_cast<int>(best.x0 + best.length - 1), candidate_end) - start_pos + 1; // 新片段长度
        int new_y0 = best.direction == 1 ? best.y0 + (start_pos - best.x0) : best.y0 - (start_pos - best.x0); // 新 y 坐标
        int correct = compute_diagonal_matches(start_pos, new_y0, new_len, best.direction, M, N); // 匹配数
        int distance = new_len - correct; // 错误数

        if (new_len >= length_thresh && distance / static_cast<float>(new_len) < rate) { // 满足长度和错误率要求
            result.push_back({static_cast<int16_t>(start_pos), static_cast<int16_t>(new_y0), 
                            static_cast<int16_t>(new_len), static_cast<int8_t>(best.direction), 
                            static_cast<int16_t>(distance)});
        }
        covered_end = candidate_end; // 更新覆盖结束位置
    }
    return result;
}

// 合并具有容差的片段
// 参数：segs（输入片段），max_gap（最大间隙），max_error_rate（最大错误率），M，N（序列长度）
// 返回：合并后的片段集合
std::vector<Segment> merge_with_tolerance_np(const std::vector<Segment>& segs, int max_gap, float max_error_rate, int M, int N) {
    if (segs.empty()) return {};

    // 按对角线分组
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
                curr.length = static_cast<int16_t>(merged_len);
                curr.distance = static_cast<int16_t>(gap + curr.distance + s.distance);
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
// 参数：matches（匹配片段集合），M，N（序列长度）
// 输出：匹配点总得分
void get_detail_np(const std::vector<Segment>& matches, int M, int N) {
    int point = 0; // 总得分
    for (const auto& seg : matches) {
        if (seg.length <= 0) continue; // 跳过无效片段
        int x0 = seg.x0, y0 = seg.y0, length = seg.length, direction = seg.direction;
        int x1 = x0 + length - 1; // 结束 x 坐标
        int affective = compute_diagonal_matches(x0, y0, length, direction, M, N); // 匹配点数
        //if (direction == 1) {
        //    int y1 = y0 + length - 1;
        //    std::cout << "(" << x0 << "," << x1 + 1 << "," << y0 << "," << y1 + 1 << ")\t";
        //} else {
        //    int y1 = y0 - (length - 1);
        //    std::cout << "(" << x0 << "," << x1 + 1 << "," << y1 << "," << y0 + 1 << ")\t";
        //}
        //std::cout << ";length = " << length << "\taffective length = " << affective << "\t";
        if (affective / static_cast<float>(length) > 0.9) { // 匹配率超过 90%
            //std::cout << "True" << std::endl;
            point += affective; // 累加匹配点
        } else {
            //std::cout << "False" << std::endl;
        }
    }
    std::cout << "final score = " << point << std::endl; // 打印总得分
}

// 主函数，寻找最佳匹配片段
// 参数：query（查询序列），reference（参考序列）
// 返回：最佳匹配片段集合
std::vector<Segment> find_best_matches_np(const std::string& query, const std::string& reference) {
    query_ptr = &query; // 设置全局查询序列指针
    reference_ptr = &reference; // 设置全局参考序列指针
    int M = query.length(), N = reference.length(); // 获取序列长度

    // 执行一系列处理步骤
    auto final_segs = find_diagonal_segments_np(M, N, 5); // 寻找初始对角线片段
    final_segs = merge_with_tolerance_np(final_segs, 3, 0.08, M, N); // 合并具有容差的片段
    final_segs = minimal_interval_cover2_np(final_segs, 0.085, 30, M, N); // 最小区间覆盖
    final_segs = fill_in_blanks_global_np(final_segs, 0.075, 15, M, N); // 全局填充空白
    final_segs = chose_segs_np(final_segs, 30); // 选择长度 >= 30 的片段
    final_segs = extend_start_backward_np(final_segs, 0.085, M, N); // 向前扩展起始
    final_segs = extend_end_backward_np(final_segs, 0.085, M, N); // 向后扩展结束

    // 按 x0 升序排序最终片段
    std::sort(final_segs.begin(), final_segs.end(), [](const Segment& a, const Segment& b) {
        return a.x0 < b.x0;
    });

    return final_segs;
}

int main() {
    std::string reference = "TATTATAGTCTTCATTCTGTGTATTAGATTACTAAAGCATATTACTTCTGTCTAAATGAAATTTTGTATCCTTTAGTCAGCATCTTCCCTTTCTCCATCCACTTTTCTCTCCCAGCCTCTGGTAATCAACATTCTACCACAATTCTGTGAGTTCTACTTTTTTAGATTCCCCATGTAAGTAAGATCATGCTGTATTTGTCTTTGTGTGCCTGGCTTATCTGAGTTAACGTAATGTCCTCCAGGTTCATTCATGTTGCAAATGATAGAATTTCCTTCTTTATCAAGACTGGATAGTATTCCATTGGTGTATATATACTATATTTACTTTATGCATTCGTCTAGACACCTAGATTGCTTCCAAATCTTAGCTATTGTGAATAGCACTGCAGTTAACATGGGAGTGCAGATATTTTTTTAGCATACTGATTTCAATCCTTTGGCCCAGAAGTGGGATTACTAGATTATATGGTAGTTCTATTTTTAGTTTTTTGAGAGACCTTCATACTGTTCTCCATAATAGCTGTGATAATTTACATTCCCACCAACAATGTACAAGTTCCCTTTTTGCTACACACTCACCAACGCTTGTTATCTTTCATCTTTTTGATAGTCTTTCTAACAGGTGTGAGGTGATCCCTCATTGTGGCTTTAATCTGCATTTCCCTGGTGAGTGGTGATGATGAGCATTTTAATATATATCTGTTGGCCATTTGTACATCTTTTGAGCAATGTTTAGGTCCTCTGCACATTTTTAAATTGGGTTCTTTGTTTTCTTGCTATTGAGTTGAGTTCTTTGTCTATTTTAGATATTAGCGCCTTATCAGATATGTAGTTTGTAGATATTTTCTCTCAATCCATGGCACATCTTTTGCTCTGTTGTTTCTTTGCTGTGTATAATCTTTTCAGTTAGATGCAATCTCATATGTTTTTGCTTTTGTTGTCAGTGCTTTTAGGGTAATATTTAAGAAATCTTTGCCCAGACCAGTGTTATGGAGCATTTTCCCTATGTTTTATTTCAGTAGCTTTACAGTTTCAGGTTTTATGTTTAAGTCATTAGTCCATTTTTAGTTGATTTTGGTGTATGGTGTGAGATAAGGGTCTAATTTCTTTCTTTTGCATGTGGTTAACCTGTTTTCCCAGCACAATTTATTGAGGATTGTCCTCTTTCCATTGTATATTCTTGGCACGTTTGTTGTAAATAATTGACCACAAATGTGTGGGTTTACTTCTGGGCTCTCTACCCTGTTTGATTGGTTAGTTGGTCGGTTTTTATGCTGTGCTGTTTTGGTTACATTAGCTTTGTAACAGGTTTTAAAATCAAGTACTGTGATACCTCTAGGTTTGTTCTTTTCGCTTTGGCCATTTGGGGTTTTTTGTGGTTCCATATGAGTTTTAGGATTGTTTTATTCTGTGAAGAATGACATTGGAATTTCGTTAGGCATTGCATTGAATCTGAATATCACTTTGGGAAGTATGAATAGTTTAACAATATTCTTGCATTCCATGAATATGGATAGCTTTCCATGTGTTTGTGTCATCTACTATATCTTTCATCAGTGTTTTGTATTTTCTAGTATATACATCTATTGCCTCCTTAAATTTACTTCTAAGTTGCTTTTTCTCGATGCTACTGTAAACGGGATTGAGTTCTTAATTTCTTTTTCAGGTACTTCATTCTTAGAATAGAGAAACACTGTAGATTTTGTATCCTGCAACTTTACTGAATTTGTTTATCAGACTTAATGGTTTTTTGGTGAAGTTGTTAGCATTTTCTATATGTGAGACCATGACATCAGCAAACAGATGATTTCTCTTCTTCCTCTGATATGGATGCCTTTATTTCTTTCTCTTGCCTAATCGCTCTTGCCAGGACATCTAACTCTGTTGAATAGCAGTGGCAAGAGCGGGCATTCTTATCTTGTTTTTGATCATAGAGGAAAAGCTTTTACCTTTTTCCAGGGAGCATGGCAGCTGTGGGCTTGTTACATATGGCCTTTATTGTGTATAAACACAATATATTCCTTCTATACCTAATGTGTTGACAGTTTTTATCATGAAATAATTTTGAGTTGTGTCACATGCTTTTTCCATATCTATTTTTTTTTATTATTATACTTTAAGTTTTAGGGTACATGTGCAAAACGTGCAGGTTTGTTACATATGTATACATGTGCCATGTCGGTGTGCTGCACCCATTAACTCGTCATTTAACATCAGGTATATCTCCTAATGCTATCCCTCCCCACCTCCCCCACCCCACGACAGGCGCCGGTGTGTGATGTTCCCCTTCCTGTGTCCATCTGTTCTCATTGTTCACTTCCCACCTGTGAGTGAGAACATGCGCTGTTTGGTTTTTTGTCCTTGTGATAGTTTGCTGAGAATGATGGTTTCCAGCTTCATCCATGTCCCTACAAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGAATATGTGCCACATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGCTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAATCGTACGTGGGCATGTGTCTTCATAGCAGCATATTTTGTAATCCTTTTGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACACTCCCACCAACAGTGTAAAAGTGTTCCTGTTTCTCCACATCCTCTCCAGTACCTGTTGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGGTATCTCACCGTGGTTTTGATTTGCATTTCTCTGATGGCCAGTGATGATGAGCATCTTTTCACGTGTCTTTTGGCGGTATAAATGTCTTCTTTTGAGAAGTGTCTGTTCATATCCTTTGCCCACTTTTTGATGGGGTTGTTTTTTTCTTGTAAATTTGTTTGAGTTCATTGTAGATTCTGGATATTAGCCCTTTGTCAGATGAGTAGATTGCAAAAATTTTCTCCCATTCTTTAGGTTGCTGTTCACTCTGATGGTAGTTTCTTTTGTTGTGCAGAAGCTCTTTAGTTTAATTAGATGGCGTTTGTCAATTTTGGCTTTTGTTGCCATTGCTTTTGGTGTTTTAGACATGAAGTCCTTGCCCATGCCTGTGTCCTGAATGGTATTGCCTAGGTTTTCTTCTAGGGTTTTTATGGTTTTAGGTCTAACATGTAAGTCTTTAATCCATCTTGAATTAATTTTTGTATAAGGTGTAAGGAAGGGATCTGGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACCATTTGTAAAATAGGGAATCCTTTCCCCATTTCTTGTTTTTGTCAGGTTTGTCAAAGATCAGATAGTTGTAGATATGCGGCATTATTTCTGAGGGCTCTGTTCTGTTCCATTGGTCTATATCTCTGTTTTGGTACCAGTACCATGGTGTTTTGGTTACTGTAGCCTTGTAGTATAGTTTGAAGTCAGGTAGCGTGATGCCTCCAGCTTTGTTCTTATGGCTTAGGATTGACTTGGCAATGCAGGCTCCTTTTTGGTTCCATATGAATTTTAAAGTAGTTTTTTCCAATTCTGTGCAGAAAGTCATTGGTAGCTTGATGGGGATGGCATTGAATCTATAAATTACCTTGGGCAGTGTGGCCATTTTCACGATATTGACTCTTCCTACCCATGAGCATGGAATGTTCTTCCATTTGTTTGTATCCTCTTTTATTTCATTGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCACGTCTCTTGTAAGTTGGAATCCTAGGTATTTTATTCTCTTTGAAGCAATTGTGAATGGGAGTTCACTCATGATTTGGCTCTCTGTCTGTTATTGGTGTATAAGAATGCTTGTGATTTTTGCACATTGATTTTGTATCCTGAGACTTTGCTGAAGTTGCCTGTCAGCTTAAGGAGATTTTGGGCTGAGGCAATGGGGTTTTCTAGATATACAATCATGTCATCTGCAAACAGGGACAATTTGACTTCCTCTTTTCCTAATTGAATACCCTTTATTTCCTTCTCCTGCCTGATTGCCCTGGCCAGAACTTCCAACACTATGTTGAATAGGAGTGGTGAGAGAGGGCATCCCTGTCTTGTGCCCGTTTTCAAAGGGAATGCTTCCAGTTTTTGCCCATTCAGTATGATATTGGCTGTGGGTTTGTCATAAATAGCTCTTATTATTTTGAGATACGTCCCATCAATACCTAATTTATTGAGAGTTTTTAGCATGAAGGGTTGTTGAATTTTGTCAAAGGCCTTTTCTGCATCTATTGAGATAATTATGTGGTTTTACTCGTTGGTTCTGTTTATATGCTGGATTACATTTATTGATTTCTGTATGTTGAACCAGCCTTGCATCCCAGGGATGAAGCCCACTTGATCATGGTGGATAAGCTTTTTGATGTGCTGCTGGATTCGGTGTGCCAGTGTTTTATTGAGGATTTTTGCATCGATGTTCATCAGGGATATTGGTCTAAAATTCTCTTTTTTGGTTGTGTCTCTGCCAGGCTTTGGTATCAGGATGATGCTGGCCTCATAAAATGCGTTAGGGAGGATTTCCTCTTTTTCTATTCATTGGAATAGTTTCAGAAGGAATGGTACCAGCTCCTCTTTGTACCTCTGGTAGAATTTGGCTGTGAATCCGTCTGGTCCTGGACTTTTTTTGTTTGGTAAGCTATTAATTATTGCCTCAATTTCAGAGCCTGTTATTGGTCCATTCAGAGATTCAGCTTCTTGCTGGTTTAGTCTTGGGAGGGTGTATGTGTCCAGGAATTTACCCGTTTCTCCTAGATTTTCTAATTTATTTGTGTAAAGGTGTTTATAGTATTCTCTGATGGTAGTTTGTATTTCTGTGGGATTGGTGGTGATATCCCCTTTATCATTTTTTATTCCGTCTATCTGATTCTTCTCTCTTTTCTTCTTTATTAGTCTTGCTAGCAGTCTATCAGTTTTGTTGATCTTTTCAAAAAGCCAGCTCCTGGATTCATTGATTTTTTTGAAGGGTTTTTTGTGTCTCTATTTCCTTCAGTTCTACTCTGATCTTATTTCTTGCCTCTGCTAGCTTTTGAATGTGTTTGCTCTTGCTTCTCCAGTTCTTTTAATTGGGATGTTAGGGTGTCAATTTTAGATCTTTCCTGCTTTCTCTTGTGGACATTTAGTGCTACAAATTTCCCTCTACACACTGCTTTGAATGTGTCCCAGAGATTCTGGTATGTTGTGTCTTTGTTCTTGTTGGTTTCAAAGAACATCTTTATTTCTGCCTTCATTTCGTTATGTACCCAGCAGTCATTCAGGAGCAGATTGTCCAGTTTCCATGTAGTTGAGCGGTTTTGAGTGAGTTTCTTAATCCTGAGTTCTAGTTTGATTGCACTGTGGTCTGAGAGACAGTTTGTTATAATTTCTGTTCTTTTACATATGCTGAGGAGTGCTTTACTTCCAAATATGTGGTCAATTTGGAATAGGTGTGGTGTGGTACTGAGAAGAATGTGTATTCTGTTGATTGGGGGTGGAGAGTTCTGTAGATGTCTATTAGGTCCGCTTGGTCCAGAGCTGAGTTCAGTGCCTGGATATCCTTGTTAACTTTCTGTCTCACTGATCTGTCTAATGTTGATAGTGGGGTGTTAAAGTCTCCCATTATTATTGTGTGGGAGTCTAAGTATCTTTGTAGGTCTCTAAGGACTTGCTTTATGAATCTGGGTGCTCCCGTATTGGGTGCATATATATTTAGGATAGTTAGCTCTTCTTGTTGAATGATCCCTTTACCATTATGAAATGGCCTTCTTTGTCTCTTCTGATCTTTGTTGGCTTAAAGTCTGTTTTATCAGAGACTAGGATTGAAACCCCTGCCTTTTTTTGTTTTCCATTTCCTTGTAGATCTTCCTTCATCCTTTATTTTGAGCCTATGTGTGTCTCTGCACGTGAGATGGGTTTCCTGAATACCGCACACTGATGGGTCTTGACTCTTTATCCAATTTGCCAGTCTGTGTCTTTTTATTGGAGCATTTAGCCCATTTACATTTAAGGTTACTACTGTTATGTGTGAATTTGATCCTGTCATTATGATGTTAGCTGGTTATTTTGCTCATTAGTTGATGCAGTTTCTTCCTAGCCTTGATGGTCTTTGCAGTTTGGCATGTTTTTGCAGTGGCTGGTACTGGTTGTTCCTTTCCATGTTCAGTGCTTCCTTCAGGAGCTCTTTTAGGGCAGGCCTGGTGGTGACCAAATCTCTCAGCATTTGCTTATCTGTAAAGGATTTTATTTCTCCTTCACTCATGAAGCTTAGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTGTCTTTAAGAATGTTGAATATTGGCCCCCACTCTCTTCTGGCTTGTAGAGTTTCTGCCGAGAGATCCGCTGTTAGTCTGATGGGCTTCCCTTTGTGGGTAACCCGACCTTTCTCTCTGGCTGCCCTTAACATTTTTTCCTTCATTTCAACTTTGGTGAATCTGAGAATTACGTGTCTTGGAGTTGCTCTTCTCGACGAGTATCTTTGTGGCGTTCTCTGTATTTCCTGAATTTGAATGTTGGCCTGCCTTGCTAGATTGGGGAAGTTCTCCTGGATAATATCCTGCAGAGTGTTTTCCAACTTGGTTCCATTCTCCCTGTCACTTTCAGGTACACCAATGAGACGTAGATTTGGTCTTTTCACATAGTCCTATATTTCTTGGAGGCTTTGTTCGTTTCTTTTAATTCTTTTTTCTCTAAACTTCTCTTCTCACTTCATTTCATTCATTTCATCTTCCATCACTGATACCCTTTCTTCCAGTTGATCGAATCGGCTACTGAGGCTTGTGCATTCGTCACGTAGTTCTCATGCTGTGGTTTTCAGCTCCATGAGGTCCTTTAAGGACTTCTCTGCATTGGTTATTCTAGTTAGCCATTTGTCTCGTTTTTTTCCAAGGTTTTTAACTTCTTTGCCATGGGTTCGAACTTCCTCCTTTAGCTCGGAGTAGTTTGATCGTCTGAAGCCTTCTTCTCTCAACTCGTCAAAGTCACTCTCCGTCCAGCTTTGTTCCGTTGCTGGTGAGGAGCTGTGTTCCTTTGGAGGAGGAGAGGCACTCTGATTTTTAGAGTTTCCAGTTTTTCTGCTCTGTTTTTTCTCCATCTTTGTGGTTTTATCTACCTTTGGGTCTTTGATGATAGTGACGTACAGATGGGGTTTTGGAGTGGATATCCTTTCTGTTAGTTTTCCTTCTAACAGTCAGGACCCTCAGTTGCAGGTCTGCTGGAGTTTGCTGGAGGTCCACTCCACACCCTGTTTCCCTGAGTATCAGCAGTGGAGGCTGCAGAACAGCGGATACCGGTGAGCAGGAAATGTTGCTGCCTGATTGTTCCTCTGGAAGTTTTGTCTCAGAGGAGTACCCGGCCATGTGAGGTGTTAGTCTGCCCCTACTAGGGGGTGCCTCCCAGTTAAGCTACTCGGGGGGGTCAGGGACCCACTTGAGGAGGCAGTCTTTCTGTTCTCAGACCTCCAGCTGCGTGCTGGGAGAACCACTACTCTCTTCGAAGCTGCAGACATTTAAGTCTACAGAGGATTCTGCTGCCTTTTGCTTGGTTATTCCCTGCCTCCAGAGGTGGAGTCTACAGAGGCAGGCAGGCCTCCTTGAACTGTGGTGGGCTCCATCCAGTTCCAGCTTCCCTGTGGCTTTGTTTACCTACTCAAGCCTTGGCAATGGTGGGCGCCCCTCCCCCAGCCTCGCTGCTGCCTTGCAGTTTGATCTCAGACTGCTGTGCTAGCAATGAGCGAGGCTCCGCAGGTGTAGGACCCTCCGAGCCAGGCGTGGGATACAATCTCCTGGTGTGCCATCTGCTAAGACCATTGGAAAAGCGCAGTATTAGGGTGGGAGTGACCCAATTTTCCAGGTGCCATCTGTCATCCCTTTCTCCATCTAGGAAAGGGAATTCCCTGACCCCTTGCGCTTCCTGGGTGAGGTGATGCCTCGCCCTGCTTCGGCTCACGCTGGGTGCACTGTACCCACTGTCCTGCACCAACTTTCCGACACTCCCCAGTGAGATGAACCCAGTACCTTAGTTGGAAATGCAGAAATCACTCGTCTTCTACGTCGCTCAAGCTGGGAGCTCTAGACTGGAGCTGTTCCCATTCAGCCATCTTGGCTCCACCAATCCATCTTACTTTTCGATGCATTTCAAAGTTGCAGATGACAGTACACATCATCCCATTATATCTCAGCATGAACATCAATAACAAGAATTCAATGTATATTTCTATTAATATGTATTTACTTTTTTCAGTAGGTTTTATATACAATAAAAGACACACCATTTGATAAGTTTTAACATATACATATAGCACTGAACTCAAACTCCAAACAAGAAATAAAATACTTCCTATCACCCAGCAAGTTCCTCTAAGACTTTTCGCAATTAATCCCCACCCCACTATCGAAGACAGCACTGCTTTTTCTTCCATTATGTATCAATTTTGCCTGTTCTAGAACTCGATGTAAATGGAATACAGAATGTACTCTTTTGTGTAAGGTTTCTTTCACCTACCATGTATTTGAAATTCATTATGTTGTTGATGTACCAAGACAGTGTTCCTTATCATTGCTAAGTTGTATTCCATTTATGACTATACCTATTTGTTGATCCAGTCTCCAACTGGTGGGTATATGAGCTGTGTCCAGTTTCGGCCTATTATAAACACTTTTGAACAAATCTTGTAGATGTGTTTTCTTTTTCTATTAAGTAAATAACTAGGAGTAGAATTGCTAAGTCATGGGTTAGGGTAGAAAAACTCTTAATTTTATAAGAGTTCTCATTGCTCTACATATATGGTAGTATTACTCTTTAAGCTAATGGGACAAGTAAAGCAGCACCTCCTTGCAGTTTACTTTTGTATTTCCTTAATAATTGACGATGTTTAGTAACATTTTCATGTGTTCACTGGCCATTCAATATAGCATCATTTATGAAGTGTGTGTACTAATCCTTCTCCAATTAATCAAGTTGTATTTTCATTGTTCAGTTGTAGAATTTCTTTACATAGTCTTAATACAAGTCCTTTACAGATACAGGTTTTGCTATTATTTTCTTCATGTCTCTTAATGAGCAGAAGTTTTTTATTTTGAAGAAGTTTATAGAGTGTTAGGTGTTTTCTTTAAAACAGTTGCTGCTTTCTGTATCCTATCAAGGAAACCTTTGCCTGTTCTTTGTTTATTCTAAAAGCTTTATACTTTAGCTTCAACATATAGATCTATGATTCATCTTGATTTATCTTTTGGGTACCATGTGAGGTAGGGAGTCAAGGTGAACTTTCATCCATACAGATGTCTAGTTGTTCCAGCACTTTTTATTTAAAAGCCTTTCTTTTCCCTCATTGGATAGCTGTCTTTGTTAAAAATAATTTATCAGTTATCTGTTCTGCTCCATTTATATGCTTTTATATCCTTATGCTGATATCACAGTCTTCATTACTGTATTATTAGAGTAAATCTTGAAATAAGGAAATTTAAGTTTTCCACCTTGTAAATCATCAGGATTGCTTTGTATAACATAGGTCCTTTCTATCTTCATATAAACGTTATAATCAGCTTGTCAATTTCTACCAAAAGAGAAAAGTCCTGCTAAGAAAATGATTTGAACCACATTAAATTCATAGATCGACTAGGGATAACTGTGCTGTTTAACTTCCAAATATTCAGATTTTTCCTGAATATCTTATTGCTATTGGTTTCTACTTTAAATTCTACTGTAGTTACAGAGTATACTCTTTAGGATTTCAAACTTTTGAAATGTAGTGAGACATGTTTTAGGACCCATATATTATTTATCTCGGTATGCATATCATGTGTACTTGAATAGAATGTACACTCTATTTTAGTTGCCCATGATATTCTACAAATGTCAGTTAGGTCAAGGTGGCTGAAGTATTTTTACAGATCTTCTATATTCTAACTGATTTTTTGTCTAGGTTTTCTATAAGTTACAGAGTAGTGTTAAAATTTCCACCTATGACTGTAGATCTTATTAGGTATATACAAATTTATGATTACTATGTATTTCTAATGAGTTTTATCACTATGAAATGCCCTTCTTTATCTCCAGTACTGTCTTGAAGTCTACTTATCTAATATTAATATAACAACTCCATCTTATACTTACTGTTTGCATGCTACACCTTTCCTCATACATTTACTTTCAATCTATCTGTGTCATTATTTTTAAAATGCATCTCTGGTAGACAGAATACAGGAAGTCTTGTTTTCTATCCAATGTATGAATCTCTACCTTCAAATAAAACACTTAATCTATTTTCAATTAAGGTAATTATTCATAGACCTGACCTAGGTTTGTAAGTTTGCTATTTGTTTTCTTTACATCCTGATATGGTTTGGCTGTGTCTCCACCCAAATCTCATCTTGAATTGTAGTTCCCATAATCCTCACACGTCGTGGAAGGGACCTGGTGGGAGGTAACTGCATCATGGAGCAGTTACCCTCATGCTATTCTCGTGATAATGAGTGAGTTCTAAGGAGATCTGATGGTTTTATAAGGGGCTTTCCCCCAACTTTGCTCTCATTCTGTTCCTTGCTGCCACCCTGTGAAGGACATGTTTGCTTCCCCTTCCTCCATGATTGTAAGTTTCCTGAGGCCTCCCCAGCCATGCTGAACTGGGTCAATTAAACCTCTTTCCTTCATAAATTACCCAGTTTCGAGTATGTCTTTATTAGCAGTATGAGAACAGACTGACACACATCTCATTTGTTTTATATTCTGCTGTTCCTTCCTTTATACCTTTTTTCATTTTCTCTAAATACATTTCAGAATTCCATTTTGACTTCTCTTGTTTTTTAGTTATCACCTCTTTGCATTTTTTTCTTTTCCAGTGTTTGTTCTAGGAATTATAATGTATATCTTTAACTTTGCAGTGTACTTAAAGCTAATTCTACCCCTTCATAAAATAAAGACCACAGTAAGCACACAGCTCCATCAACCTAGACTCATCCTTTATACACAGCAGTCATATATAAATGCCATATTGTAACATTATAATTTCTGCTTTAAACAGTCACAACTGTATATTTCTTAATGGAGTGGGGAAGACTTACATTTTTACCAATTGCAGTGCTCTTCATTCCTGAAGTTCTGAGTTTCAGCTGGTAACACTTTCCTTCAGCTTGAAGTATCTCCTTAAACATTTCTTGTAGTATAAAGATTTTCTTAGCTAGAAAGTATGTTTTCTTATTTGGCTACAGAATTCTGAGTTGATTTTTATTTTCCCAGTACTTTAAAGGTGTTTATAGCCTGCTGGCCTGCATTATTTCTGATAAGAAGTCAGGTCTGTCAAAGATCAGATGGTTGTAGATATGTGGCATTATTTCTGAGGGCTCTGTTCTGTTCCACTGATCTATATCTCTGTTTTGGTACCAGTACCATGCAGTTTTGGTTACTGTAGCCTTGTAGTATAGTTTGAAGTCAGGTAGCGTGATGCCTCCGGCTTTGTTCTTTTGGCTTAGGATTGACTCGGCAATGCGGGCTCTTTTTTGGTTCCATATGAACTTTAAAGTAGTTTTTTCCAATTCTGTGAAGAGAGTCATTGGTAGCTTGATGGGGATGGCATTGAATCTATAAATTACCTTGGGCAGTATGGCCATTTTCACGACATTGATTCTTTCTACCCATGAGCACGGAATGTTCTTCCATTTGTTTGTATCCTCTTTTATTTCATTGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCACGTCCCTTGTAAGTTTGGATTCCTAGGTATTTTATTCTCTTTGAAGCAATTGTGAATGGGAGTTCACTCATGATTTGGCTCTCTGTTTGTCTGTTATTGGTCTAGAAGAATGCTTGTGATTTTTGCAAATTGATTTTGTATCCTGAGACTTTGCTGAAGTTGCTTATCAGCTTAAGGAGATTTTGGGCTGAGATGATGGTGTTTTCTAGATATGCAATCATGTCATCTGCAAACAGAGACAATTTGACTTCCTCTTTTCCTAATTGAATACCCTTTATTTCCTTCTCCTGCCTGACTGCCCTGGCCAGAACTTCCAACACTATGTTGAATAGGAGTGGTGAGAGAGGGCATCCCTGTCTTGTGCCAGTTTTCGGGAAAGGATTCCCTATTTAATAAATGGTGCTGGGAAAACTGGCTAGCCATATGTAGAAAGCTGAAACTGGATCCCTTCCTTACACCTTATACACAAATTAATTCAAGATGGATTAAAGACTTAAATGTTAGACCTAAAACCATAAAAACCCTAGAAGAAAACCTAGGCAATACCATTCAGGACATAGGCATGGGGAAGGACTTCATGTCTACAACACCAAAAGCAATGGCAACAAAAGCCAAAATTGACAAATGGGATCTAATTAAACTAAAGAGCTTCTGCACAGCACAAGAAACTACCATCAGAGTGAACAGGCAACCTACAAAATGGGAGAAAATTTTCGCAACCTACCTCATCTGACAAAGGGCTAATATCCAGAATCTACAATGAACTCAAACAAATTTACAAGAAAAAAACAAACAACCCCATCAAAAAGTGGGTGAAGGATATGAACAGACACTTCTCAAAAGAAGACATTTATGCAGCCAAAAAACACATGAAAAAATGCTCACCATCACTGGCCATCAGAGAAATGCAAATCAAAACCATAATGAGATACCATCTCACACCAGTTAGAATGGTGATCATTAAAAAGTCAGGAAACAACAGGTGCTGGAGAGGATGTGGAGAAACAGGAACACTTTTACACTGTTGGTGGGAGTGTAAACTAGTTCAACCATTGTGGAAGTCAGTGTGGCGATTCCTCAGGGATCTAGAACTAGAAATACCATTTGACCCAGCCATCCCATTACTGGGTATATACCCAAAGGACTATAAATCATGCTGCTATAAAGACACACGCACACCTATGTTTACTGCGGCACTATTCACAATAGCAAAGACTTGGAACCAACCCAAATGTCCAACAATGATAGACTGGATTAAGAAAATGTGGCACATATCCACCATGGAATACTATGCAGCCATAAAAAATGAAGAGTTCATGTCCTTTGTAGGGACATGGATCAAACTGGAAACCATCATTCTCAGCAAACTATCGCAAGGACAAAAAACCAAACACCGCATGTTCTCACTCATAGGTGGGAACTGAACAATGAGAACACTTGGACACAGGAAGGGGAACATCACACTCCGGGGACTGTTGTGGGATGGGGGAAGTGGGGAGGGATAGCATTAGGAGGTATACCTAATGTTAAATGATGAGTTAATGGGTGCAGCACACCAACATGGCACATGTATACATATGTAACAAACCCGCACATTGTGCACATGTACCCTAAAACTTAAAGTATAATAATAATATAAAAAAAAAAAGAAGTCAGCTACCATTCATTAGCATTTCTCAGGGTGCTTTCAAGATTTTTTTCTTTTTGTAGAGATGAAGTCTCACTTCGTTTCCCAGGCTGGTCTCAAACTCCTTAGCTCAAGCAATACCCCTCTTTGGCCTCCCCAAGTGTTGGCATTACAGGCATAAGCCACTGCAACTAGGCAGGATTTTTTGTCTATTATTTTCAACAGCTTCACTATAATGTACTTAAGTGTTATTTTATTGTATTTATCTTGTTTGGAATTCACTGAACTTCTTCAACATGTAAATTCATGTTTTTTACCAATTTAGGGAAATTTCCAGCCATAATTTCTCCAGAATTTTTTGACCCAATTTCCCCTGGCCTCTATTTCTAGGACTCCAAGTGTTAGACCTTTTGATATTGTCCTACCAGCCACTGAGGCTAGCTTCATTTTTTGAAATCTCTTTTTCACTTTTTAAGAAAGAATAACTTCTATGAATCCAACAACTTCAATGACTCCTTTTCTTTTCCTAATCTTCTGTAATCTCTAATCAGTTTTCATTTGAAATATTGTATTTTTTACCTATAGAATTTTCATTTAGTCCTTTTACCTAGTTTCCATTTTTCAGCTGAGGTTTCTTATCTGTTTATTCATGGTGTACATATTTTCTTTATATAACCTTAAGAAGAGTTATAACAGCTGCTAATTCTAACATAGGAGTCATATCAAGGCCAGTCTCAATTAATTGGCTTTTCACTTAAGTATGGGTCACATTTTTCTGTTTGTGTATGTGTGTGTCTACTGAATATGGATTGTATTCTGGACACTGTGCATAGTACATTTTTAGAAACTCTGGATTTGATTATGTTCTTCTGAAGACTATTTAGTTTTTTTTTGCTTAATAGGTAGGTAAGTTGGCTATACACCTCTGAACTCTTTCTTGTGCAAATGAAATATCTTTAGTTATTTTAGCCTAAGACTTTTATGCAGAATTTGGGACACGCTTCTTTATGAAAGACTACTTTATACTCACTAACTTTCTAGCTGCTAGAGTCACCCCTAACTCTGTCTTGTGGTTTTTATTCAAGCTTCCGTCACCTACACAACCAATGACCGAGTTCTGTCCTCTGGGAAACATCCACAAAAAATGAAATGCTCATCTTATACCATTCTATCTTCCAAATGTCAACATCCTTCAGTTTCTTCCTACTTTCTGTCCTCCTACTGTGCTTTCACATACTTTTTTTATATTTTGTCCAAAAGTTGTCATTGTTATTTTCAGGAAAGGTAGTCCGACATGTGCCCCTCTGTAATTATTAAAAGCACAACTCCAAAGAATTGTTTATAATACCAGAAGACTACAAAAGTGACATCCATAGTGACTTTTTTTTTTTTTTTTTTTTGAGATGGAGTCTCGTTCTGTCACCCAGGCTGGAGTGCAGCGGCGCGATCTCAGCTCACTGCGACCTCCACCTCCTGGGTTCACACCATTCTTCTGCCTCAGCCTCCCGAGTAGCTGGGACTACAGGCGCCCACCAACACGCCCGGCTAATTTTTTGTATTTTTAGTAGAGACGGGGTTTCACTGTGTTAGCCAGGATGGTCTCAAATCTCCTGACCTCGTGATCTGCCTGCCTCGGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACAGCGCTGGCCGACTTATTACTTTGTTTTTTGAGACAGAGTCTGTCGCCCAGGCTGGAGTGCAGTGGTGCAATCCCGGCTCACAGCAAGCTCTGCCTCCCGGGTTCACGCCATTCTCCTGCCTCAGCCTCCTGAGTAGCTGGGACTACAGGCGCCCGCCACCATGTCCGGCTAATTTTTTTTGGTATTCTTTAGAGACAGGGTTTCACTGTGTTAGCCAGAATGGTCTTGACCTCCTGACCTCGTGATCCGCCCACTTCAGCCTCCCAAAGTGCTGGGATTACAGGCGTGAGCCACTGTGCCTGGCCAACTTACTACTTTTTAACATTACACAGTTGTAACAATATTTCCTAAAGCAATCATTTTCTGTTATAGTCCTTCTTATATCAGAAAACAGAAGTATTTGGGAAATGAAAAGATGGAGAAGAAAGAGGGGAATGAAGGGAAATAAGATGACCCAGTGATACTGAAAGGGGAAGAAGGTTCATGCAAAGAAAAGGAACAATCTGACAAACTCCTTGTCCTGGCATCTTTTAAGCTAAAGTTCAGCTATACAAGAAGGGAGATGATGAAAAGATAAAATATAATTTCTGGCACAGATGCAACAATTGATATTGGTACCTACAGTTACTCACTTTCATAAAGAGAATAAAAGCCCGAAAATGCTATACTAAGATAATGAACAGATACAGTATATAATACATGAGTTATAATATAAAATAACAGTAAAAAAATACAGACTCCTCAATAAAATCTCCAACCTTTTCTTTTTACCATAAGAACACACTGGAATATAGACTTCAGTACAATGTTTCTATAATTAAATCTTTAAAATATAATCATACAAAAAGCTTGACTAAGAACTCCAGATTTCCCTTTACAAGCACTAAATTTTATCAGTTTTTATATAAACATATATGAGAGGTTTCTAAAATGTATGATGATCAGAATTAAACTCTTTCAACTATAGGTAACTTATACGCAATCAGAGGATTCTGGGAAGATGGTGGAGTAGGAAGCACCACGCATGTGTTGCTCCACCTAGAAGATCATTGCACTGGCAGAATCTGTTCTATTTAACTATTTTGAAACTCTGGAATCTACTGAAGGCTTGTTAATGTCCAGGGGAAGGCCTGAATATAATTTTAATCACTTGCAGGTATTAGCACTAGCACAGTAGCAGCTATCCATTTGCCACCCTCAGCCCCAAGGCAGGCAGTTGTGCATGCATTCCAAGGAGCAGCTTGCGCGCAGCTTACAGGAACCAAGGACAGGCAAAAAGGACCCTCTCTTCCAAATATAGGGCATGTGTGCCCTGATCACTGAATGCTACTTCTGATCACAGAAATGGGACAAAGAGGTGGACACACATTGTTGCTTCATCTCCACTCACTGTTACAAGCCCTCTCCCACCCCCATCCCCACCCCCACCCACAAAGTAACTGCTGGGGATTTAAAGGATTTAGTTTTTCTCTTTTTCCCTTTTGGGAGTCAGATGTTAAAGACTACGGCATTTGAAAGCAACTGTATATAAGGGAAAAATTAGAAACTTCCTGCACATTCCTGGAGAAAAACTCAGGCTTAGAAAAGACTCAGAGAATCTTAAGCTTACACCTCAGGCTGATCCTTCACATACAGACAAGCCTCAAACAATAAAAATAAACAAAAACAGTAACATTAGAAAAAAAAAAAAAAAACAGCAGGCCGGGCATGGTAGCTCATGCCTATAATCTCAGCACTTTGGGAGGCTGAGGTGGGCAGATCACCTGAGGTCAGGAGTTTGAGACGAGCCTGGGCAACATGGCGAAACCCCATCTCTACTAAGAAGCACAAAAATTAGCTGGGCATGGTGGCACACACCTGTAATCCCAGCTACTTGGGAGGCTGAGGCAGGAGAATCGGCTGAACCCAGGAGGTAGAGGTTGCAGTGAGCCGAGATCACACCACTGCACTCCAGCCTGGGCAACAGAGTGAGACTCTGCCTCAAAACAAAACAAAACAAACAGCAAATGCTGGGGAAGAGGGAGAATCTAATTTCCAGACTTACCACATCATTATATTCAAATGCCCTGTTTTCAACAAAAAAATAACAAGATATACAGGAAAGTATGACACATTCAAAAGAAAAAAAAATAAAACAACACAACTGTCCTTGAAAAAGACCTGATGGCAGACCTTCCAGAGAAAGACTTTATGACAACTGTTTTAAAGTTGTTCAAAGTACTAGGGCGGGCACAGTGGCTCATGACTGTAATCCCAGCACTATGGAAGGCCAAGGCGGGAAGATCACCTGAGGTTGGCAGTTTGAGACCAGCCTGGCCAACATGGCAAAACCCCATCTCTACTAAAAATACAAAAATTAGCCGGGCATGGTGGCACACGCCTGTAATCTCAGCTACTAGGGAGGCTGAGGCATGAGAATCGCTTGAACCCGGGAGGCAGAGGTTGTAGTGAGTTGAGACTGTGCCACTGCACTCCAGCTCGGGTGACAGAGTGAGACTCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAACGAAATGAGAAAAGAGGCTGGGTGCGGTGGCTCACACCTGTAATCCCAGCATTTTGGGAGGCCAAGGTGGGTGGATCACCTGAGGTCAGGAGTTCGAGTCCAGCTGGGCCAACATTGTGGAACCCCATCTCTACTAAACATACAAAAATGAGCCAGGCGTGGTGGCACATGCCTGTAATCCCAGCCTCAGGAGGCTGAGGCAGGATAACTGCTTGAACCCAGGAGGCAGAGGTTGCAGTGAGCCAAGATCATGCCACTGCACTCCAGCCTGGGCGATAGAGACTCCATTAAAAAAAAAAAAGAAGAAAGAAAGAAGAAAGGAATTTAAACATTTCATGAAAAAAATCAAGTAAACACAAAAGAAGACAGTAATGCAGGAAATGAGGGCCCCCACAAAAAGCAGTAAGTTATGTAAACAACAAACAACAAAATGACAGAAATCCTTCTGTACAAATAATTTAAATTCAAATGGATTAAACTCTCCAATCAAAAGACAGAGATTGACAGAATGGAGGAGGATTCAACCATATGTTGCCTATACTAGACTTACTTTAGACTGACAGACACAGTTTGAAAATTAAAAGACAGAAAAAGATATTCCATGCAAATACTAACCAGAACAGAGTAGAGGTGGCTATACTAACAACAGAGAAAACAGTCTAAATCAAAAAAGGTTACAAGAGACAAAGGACATTATATATGAATACATGATTTGATACAACAAGAAGACAAAAATTATAAACATTTATGCACCAATTAAAAGACCTCAAAAATACATGAAGGAAAAACAGACAACTGAAGGGAGAAACAGACCTTTCTAAATAACAGATGATTAAATCAATACCCTACTTTCAATAGTGGATAGAATAACCAGATCGAAGATAAGTAAGGAAACAGAGGACCTCAGACAATAACCAACTAGATCTGACAGACATATACACAACACTGTACCCAACAATAGCATACATATTCTTCTCAAGTGAATATAGAACATTTTCCAGGATAGACCATATGTTAGGTCACAAATTAAGTAAAAAATGTGGTATAAACATACAATGGAATATTATTCAGCCTTTAAAAAGAAGGAAATTCTAGCCGGGCTCGGTGGCTCATGCCTGTAACCCCAGCACTTTGTGAGGCTGAGGTGGGCGGATCACTTGAGGTCAAGAGTTCAAGACCAGCCTGGCCAACATAGTGAAGCCCTGTCTCTAGTAAAAATACAAAAATTAGCTAGGCGTGGTGGCGTGCGCCTATAATCCCAGATACTAGGAGGCCAAGGCACAACAATCACTTGAACCCAGGAGGCGGAGGTTGCAGTGAGTCAAGATCGCACCACTGCACTCCAGCCTGGCCAACAGAGTGAGACTCCATCTCAAAAAATAAAAATAAAAAATTCTAACAACATGGATTAACCTTAAGGATATTATGCTAAGTGAAATAAACCAGTCACAGAAAAACAAATACTGTATGATTCTATTCATATGAGGTAATTACAGTAGTCAAAATCATACACAAAGTAGAATGGTGGTTGCCAGGAGCTAGGGAAAGACAGAAGGAGGAGTTACTGTTTAATGAGTACAGAATTTCAATTTTACAAGATGAAAAGAGTTCTGGAGATAGAAGGTGGTGATGGTAGTTACACAACATTATGAAGATACTCAATACCACTGAATTGCACATTTTAAAATGGTCAAAATAGTAAATTTTATGTTATGTGTATTTTAAAATGTAAAACATTGGGAAAAAAGATGAAAACAATGGGCAGTTCTTAGGCATTAACAATTAGTACACATCTAGTCATATCATTTTATTCTAACTCACTATTAAAACCGATTAAAAGCAAAACTGATTTTCAAAAATTTATATAAATAAGTAAAAATACAGAAAATTATACATATGCACATTAAATTTAGAAATGTTTTTAAAGGTTTGTCAGATCACTTTCTGCAAGATTATAATACACTCCCATGTAGTACAAAATTTTGTTTATTCCTGTTTTTTTAAACTTTTCCGTAATTCCATCTTCTATGAGTTACTGAATCATTTTTTCTACCCTTACTGATAAGAATACAAACAAAATATCTGTTTCATAAAAATGAATCCGTATTTCAAATTTGAAAGCCACATTTACTTTTCCTTTAAAACTACATGAACTACTTTCAGATCTCTTATCTGTGCCAGCCAGCTTGTTGTTTACCTTCTCACAAATGTGCCTGAGTCTCACTCACATTACTTTACACTGATTATGTCTATTAGCCTGTATTTGCATAGGAGTCACTACATAAAGTATATAAGAATGACATAACAAAGTTCACCTAAATCCAATTTCCTTCACACTTCCTGAAAAATTAGGTTATGAGAATTCACTGCTATGTCTCTGCTGCCCTCTGTAGGTATGGTGTCTTAGTCCACTTTGCTTTGCTTTAACAGAATACCACAAGATAGGTAACTTACAAAGAAAAGAAGTTTATTTCTCATAGTTATGGAGGCTGAAAATTCCAATATCAAGGTGCTGGCTCCTGGTGAGGACCTTTGTCCTACACCTTCCCATGGCAGAAAGCAGAACGACAAGTGAGGGTGAGAGCAAGGGACAGAAAATTTTTATAACAACACATTCTTGTGGTAACAAACCCACTACTGTGATAATGATATTAATCTATTCATGAGGGCAAGGCCCTCATGGCCTAATCCTAATACAGGAACCATTAAAAGGTTCTTTTAATATATACCATTCTTTTAATGTATACCATTCTCACATTGGTATAAAGATATACCTGAGACTGGGTAATTTATAAAAGACGAGGTTTAACTGGCTCACAGTTCTGCAGGCTGTACAGGAAGAATGGTGGCATATGCTTCTGGGGACACCTCAGGGAGCTTTCACTCATGGTGGGAGGCAAAGTGGGTGCAGGCACTTCACATGGTGAAAGCAGGGGGAGGAGGAAGCAGACACTTCACATGGTGAAAGCAAGAGAGCAAAGTGGGCAGTACTACACACTTTTACGTATTTATTTATTTTAATTTTTTTTAATAGAGACGAGGTCTTGCTATGTTGGCCAGGTCGGTCTCAAACTCCTGGGCTCAACTGATCCTCCCATCTTGGCCTCCCAAAGTGCTGGGGTTACAGGCAGGAGACACAGTGCCTGACCACTACACACTTTTAAATTATCAGCTCACATAAGAACTCACTCACTATTGTGAGGACAGTACCAAGGGTTTTGGTGCTAAACCATTCATGAGAAATCCACTCCCATAATCCAGTCACCTTCCACCAGGCCCCATTTCCAACAATGGAAATTACAATTCAACATGAGATTTGGGCGGGGACACAGACTCAAACTATATTAGTCCCCTATTTTAATACTGTCACCATGGTAATTAAATTTCAACAAGAGTTTTGGAGAGGACATTCAAACCATAGCAAATGCTTATTGGAAGATACCAATAAAATAAATTATAAAAACCCAGACAAGAATAAAAAGTACTATCAACTTAATATGGAAAGTCTTATCAGTAATCAAATACATCTGTTAAATGTCCTAAAAGAAAATACCGCAGGCCAGGTGCAGTGGCTCACACCTGTAATCTCAGCACTTTGGGAGGCCGAGGTGGGTGGATCACCTGAGGTCAGGAGTTCAAGACCAGCCTGGCCAATACGGTGAAACCCTGTCTCTACTAAAAATACAAAAATTAGTTGGGCATGGTGACGGGCGCCTGTAATCCAAGCTACCATGTCATCTTTGATGAACAAGGATGCAAAAATCACCAATATACTAATAAACTGAATTCAACAGCACATTAAAAGGATCATTCACATGATCAAGTGGGATTTACCCTTGGGATGCAAAGATGGTTCAACATATGCAAATAAATAAATGTGATACACCACGTTAATAGAATGGAAGACAAAATCCATATGATCACCTGAATAGATATGGAAAGGGGTGGAGCCAAGATGGCCAAACAGGAACAGCTCCAGTCTACAGCTCCCAGCGTGAGTGACGCAGAAGACGGGTGATTTCTGCATTTCCAACTGAGGTATTGGGTTCATCTCACTGGGGAGTGTTGGACAGAGGGTGCAGCACACTGGGTGCAGTGCACCCAGCGTAAGCCGAAGCAGGGCGAGGCATCACCTCACCCGGGAAGTGCAAGGGGTCATGGAATTCCCTTTCCTAGTCACAGAAAGGAGTGACAGACAGCACCTGGAAAATCGGGTCACTCCAACCCTAATACTGCGCTTTTCTAACGGTCTTAGCAAATGGCACACCAGGAGATTATATCCCACAACTGGTTCGGAGAGTCCTACGCCCACGGAGCCTCGCTCATTGCTAGCACAGCAGTCTGAGATCAAACTGCAAGGCAGCAGCCAGGCTGCGGGAGGGGCAGCCACCATTGCCAAGGCTTGAGTAGGTAAACAAAGCCACAGGGAAGCTGGAACTGGGTGGAGCGCACCACAGCTCAAGGAGGCCTGCCTGCCTCTGTAGACTCACCTCTGGGGGCAGGGCATAACCAAACAAAACGCAGCAGAAATCTCTGCAGAATTAAATGTACCTGTCTGACAGCTTTGAAGAGAGTAGTGGTTCTCCCAGCATGCAGCTGGAGATCTGAGAACGGACACACTGCCTCCTCAAGTGGGTCCCTGACCCTCGAGTAGCCTAACTGGGAGGGACCCCCCAGTAGGGGCAGACTGACACCTCACATGGCCAGGTACTCCTCTGAGACAAAAACTTCCAGAGGAATGATCAGGCAGCAACATTTGGTATAAACCAATATTGGCTGTTCTGCAGCCTCCGCTGCTGATACCCAGGCAAACAAGGTCTGGAGTGGACCTCCACCAAACTCCAACAGACCTGCAGCTGAGGGTCCTGACTGTTAGAAGGAAAACTAACAAACACAAAGGACATGCACACCAAAACCCCATCTGTACGTCACCATCATCAAAGACCAAAGGTACATAAAACTACAAAGATGGGGAAAAAACAGAGCAGAAAAACTGGAAACTCTAAAAATCAGAGTGCCTCTCCTCCTCCAAAGGAACACAGCTCCTCACCAACAACGGAACAAAGCTGGACGGAGAATGACTTTGATGAGTTGAGAGAAGAAGGCTTCAGACGATCAAACTACTCCGAGCTAAAGGAGGAACTTCGAACCCATGGCAAAGAAGTTAAAAATCTTGAAAAAAATTAGACGAATGGCTAACTAGAATAACCAATGCAGAGAAGTCCTTAAAGGACCTCATGGAGCTGAAAACCACAGCACGAGAACTACGTGACGAATGCACAAGCCTCAGTAGCCGATTCGATCAACTGGAAGAAAGGGTATCAGTGATGGAAGATGAAATGAATGAACTGAAGTGAAAAGAGAAGTTTAGAGAAAAAAGAATTAAAAGAAACGAACAAAGCCTCCAAGAAATATAGGACTATGTGAAAAGACCAAATCTACGTCTCATTGGTGTACCTGAAAGTGACAGGGAGAATGGAACCAAGTTGGAAAACACTCTGCAGGATATTATCCAGGAGAACTTCCCCAATCTAGCAAGGCAGGCCAACATTCAAATTCAGGAAATACAGAGAACGCCACAAAGATACTCGTCGAGAAGAGCAACTCCAAGACACGTAATTCTCAGATTCACCAAAGTTGAAATGAAGGAAAAAATGTTAAGGGCAGCCAGAGAGAAAGGTCGGGTTACCCACAAAGGGAAGCCCATCAGACTAACAGCGGATCTCTCGGCAGAAACTCTACAAGCCAGAAGAGAGTGGGGGCCAATATTCAACATTCTTAAAGACAAGAATTTTCAACCCAGAATTTCATATCCAGCCAAACTAAGCTTCATGAGTGAAGGAGAAATAAAATCCTTTACAGATAAGCAAATGCTGAGAGATTTGGTCACCACCAGGCCTGCCCTAAAAGAGCTCCTGAAGGAAGCACTGAACATGGAAAGGAACAACCAGTACCAGCCACTGCAAAAACATGCCAAACTGCAAAGACCATCAAGGCTAGGAAGAAACTGCATCAACTAATGAGCAAAATAACCAGCTAACATCATAATGACAGGATCAAATTCACACATAACAATATTAACCTTAAATGTAAATGGGCTAAATGCTCCAATAAAAAGACACAGACTGGCAAATTGGATAAAGAGTCAAGACCCATCAGTGTGCGGTATTCAGGAAACCCATCTCATGTGCAGAAGACACACACAGGCTCAAAATAAAGGATGGAGGAAGATCTACCAAGGAAATGGAAAACAAAAAAAGGCAGGGGTTTCAATCCTAGTCTCTGATAAAACAGACTTTAAACCAACAAAGATCAAAAGAGACAAAGCAGGCCATTATGTAATGGTAAAGGGATCAATTCAACAAGAAGAGCTAGCTAGCCTAAATATATATGCACCCAATACGGGAGCACCCAGATTCATAAAGCAAGTCCTTAGAGACCTACAAAGATACTTAGACTCCCACACAATAATAATGGGAGACTTTAACACCCCACTATCAACATTAGACAGATCAGTGAGACAGAAAGTTAACAAGGATATCCAGGTATTGTACTCAACTCTGCACCAAGCGGACCTAATAGACATCTACAGAACCCTCCACCCCAAATGAACAGAATATACATTCTTCTCAGCACCACACCGCACTTATTCCAAAATTGACCACATAGTTGGAAGTAAAGCACTCCTCAGCAAATGTTAAAAGAACAGAAATTATAACAAACTGTCTCTCAGACCACAGTGCAATCAAACTAGAACTCAGGATTAAGAAACTCACTCAAAACTGCTCAACTACATGGAAACTGAACAACCTGCTCCTGAATGACTACTGGGTACATAACGAAATGAAAGCAGAAATAAAGATGTTCTTGGAAACCAACGAGAGCAAAGACACAACATAACCAGAATCTCTGGGACACATTCAAAGCAGTGTGTAGAGGGAAATTTATAGCACTAAATGCCCACAAGAGAAAGCAGGAAAGATCTAAAATTGACACCCTAACATCACAATTAAAAGAACTAGAGAAGCAAGAGCAAACACATTCAAAAGCTAGCAGAAGGCAAGAAATAACTAAGATCAGAGCAGAACTGAAGGAAATAGAGACATAAAAAACCCTTCAAAAAATCAGTGAATCCAGGAGCTGGTCTTTTGAAAAGATCAACAAACTTGATAGGCTGCTAGAAAGATTACTAAAAGAGAAAAGAGAGAAGAATCAAATAGATGCAATAAAAAATGATAAAGGAGATATCACCACTAATCCCACAGAAATACAAACTACCATCAGAGAATACTATAAACACCTCTACGCAAATAAACTAGAAAATCTAGAAGAAATGGATACATTCCTGGACACATACACCCTCCCAAGACTAAACCAGAAAGAAGTTGAATCTCTGAACGGACCAATAACAGGCTCTGAAATTGAGGCGATAATAATAGCTTACCAAACAAAAAAAGTCCAGGACTGGATGAATTTACAGCCGAATTCTACCAGAGGTACAAGGAGGAGCTGGTACCATTCCTTCTGAAACTATTCCAATCAATAGAAAAAGAGGGAATCCTCCCTAACTCATTTTATGAGGCCAGCATCAGCCTGATACCAAAGCCTGGCAGAGACACAGCAAAAAAAGAGAATTTTAGACCAATATCCCTGATGAACATGGATGCAAAAATCCTCAATAAAATACTGGCAAACCAAATCCAGCATCACATGGAAAAGCTTATCCACCATGATCAAGTGGGCTTCATCCCTGGGATGCAAGGCAGGTTCAACATATGCAAATCAATAAACATAATCCATCATATAAACAGAACCAACGAGTAAAACCACATGATTATCTCAATAGATGCAGAAAAGGCCTTTGACAAAATTTAACAACCCTTCATGCTAAAAACTCTCAATAAATTAGGTATTGATGGGACGTATCTCAAAATAGTAAGAGATATCTATGACAAACCCACAGCCAATATCATAATGGGCAAAAACTGGAAGCATTCCCTTTGAAAACTGGCACAAGACAGGGATGCCCTCTCTCACCGCTCCTATTCAACATGGTGTTGGAAGTTCTGGCCAGGGCAATCAGGCAGGAGAAGGAAATAAAGGGTATTCAATTAGGAAAAGAGGAAGTCAAATTGTCCCTGTTTGCAGATGATACGATTGTATATTTATAAAACCCCATCGTCTCAGCCCAAAATCTCCTTAAGCTGATAGGCAACTTCAGCAAAGTCTCAGGATACAAAATCAATGTGCAAAAATCACAAGCATTCTTATACACCAATAACAGAGAAAGAGAGAGCCAAATCATGAGTGAACTCCCATTCACAATTGCTTCAAAGAGAATAAAACACCTAGGAATCCAACTTATAAGGGATTTGAAGGACCTCTTCAAGGAGAACTACAAACCACTGCTCAATGAAATAAAAGAGGATACAAACAAACAGAAGAACATTCCATGCTCATGGACAGGAAGAATCAGTATCATGAAAATGGCCATACTGCCCAAGGTAATTTATAGATTCAGTGCGATCCCCATCAAGCTACCAATGACTTTCTGCACAGAATTGGAAAAAACTACTTTAAAGTTCATATGGAACCAAAAAAGAGCCGGTATTGCCAAGTCAATCCTAAGGTGAAAGAACAAAGCTGGAGGCATCGTGCTACCTGACTTCAAACTATACTACAAGGCTGCAGTAACCAAAACTGCATGGTACTGGTACCGAAACAGAGATATAGAGCAATGGAACAGAACAGAGCCCTCAGAAATAATGCCACATATCTACAACCATCTGATCTTTGACAAACCTGACAAAAACAAGAAATGGGAAAACGATTCCCTATTTAATAAATGGTGCTGGGAAAAGTGGCTAGCCATATGTAGAAAAGTGAAACTGGATCCCTTCCTTACACCTTATACTAAAATTAATTCAAGATGGATTAAAGACTTAAATGTTAGACCTGAAACCATAAAAACCCAAGAAGAAAACCTAGGCAATACCATTCAGGACATGGGCATGGGCAAGGACTTCGTGTCTAAAACACCAAAAGCAATGGCAACGAAAGCCAAAATTGACAACTAGGATCTAATTAAACGAAAGAGCTTCTGCACAGCAAAAGAAACTACCATCTGAGTGAACAGCAACCTACAGAATGGGAGAAAATTTTTACAGTCTGCTCATCTGACAAAGGGCTAATATCCAGAATCTACAATGAGCTCCAACAAATTTTCAAGAAAAAAACGCCATCAGAAAGTGGTCGAAGGATATGAGCAGACACTTCTCAAAAGAAGACATTTATGCAGCCAAAAGACACTTGAAAAGATGCTCATCATCCCTGGCCATCAGAGAAATGCAAATCGAAACCACAATGAGATACCATCTCACACCAGTTAGAATGGCAGTCATTAAGAAGTCAGGAAACAACAGGTGCTGGAGAGGATGTGGAGAAATAGGAACACTTTTACACCGTTGGTGGGACTGTAAACTAGTTCAACGATTGTGGAAGTCAGTGTGGCGATTCCTCAGGGATCTAGAACTAGAAATACCATTTGACCCAGCCATCCCATTACTGGGTATATACCCAAAGGATTATAAATCATGCTGCTATAAAGACACATGCACACGTATGTTTATTGGGGGACTATTCACAATCGCAGAGACTTGGAGCCAAGCCAAATGTCCAACAATGATAGACTGGATTAAGAAAATGTGGCACATATACACCATGGAATACTATGCAGCCATAAAAAAGGATGAGTTCATGTCCTTTGTAGGGACATGGATGAAGGTGGAAACCATCATTCTCAGCAAACTATCAGAAGGACAAAAAACCAAACAGCGCATGTTCTCACTCATAGGTGGGAATTGAAAAGTGAGAACACAGGGGCACAAGAAGGGGAACATCATACGCCGGGGCCTGTTGTCGGGTGGGGGCAGGGTGGAGGGATAGCATTAGGAGGTATACCTAATGTTAAATGACGAGTTACTGTGTGCAGCAGACCAACATGGCACATGTATACATATGTAACTAACCTGCACATTGTGCGCATGTACCCGAAAACGTAAAGTATAATAAAAAAAGAAGAAGAAGAGAAAAAAAAAAGAAAAGTAAGGTGGATTAATGGATAGAGAAATGTGATAAACTGTCATAGCTAAATAATAACAATTGTAAAGCCTTGGTGATGTATATATTGGTTTTTGCTAGTATCTCTAAAAAATTTCTGGTCGAAATTTGCATAATATTCTGTGTGAAAAATCTTCTCAGGGAACATAGACCTGTGTTATTGTATTTATGTAAGTAGTTGTCATATAATTCTTTTATCTTAGAGAATCTGAAGTTTCCAATGTGAGTGTACAGTGTCTCACCACCACCTGAGTTTATTTCTGGTTGTTGTAATGCAGTACAATTGATAATGAAATTTTTTTGATGTTTTTGCCGTAGACATGGAATGTATAGTGCTCTAGCTTTTACTTTTTGAGAGACAGAAAAAGGCCCTGAGCATTGGTGTTGTGTCAGTCATCCTAGGCTCAGATCTCCATTCTACCTCTCACTGAATGCTTGGCCTTAAGCAGGTATTTAACTTTTTTTTAGATCCAATTGGTTTGTCTGTTAAATGCCTTCAGTAACGTATTCCTCACATAGTTGCTTTGAAGATTGGTGAGTCATGGCACAGCACACAGCTGAACATAGAAGACATTCAAGCAATGCCAGATTCCTTTAGGGTGCAGGTAGACATTCCTAGGGTTTGTAATAGATTTAGAAGTAGAATTAGTCTTCTACTTCTTCAGATCTCCTTTGAGGTCTGAGTGACTACAGAACCCATGGACACCTTTATTTACGAAATGCTGTGTTTCTGGTGTCTAAGGCCTGATTGCATCTTGATTGGAAACTATTAGTTAAACATTTATAAGAAATAATACTGTAACTTTTGTTTAAGCTAAAAATTATTTTCTTGAATTTTAGTTTATTGGTCTTAGTTAATACTATAACCTTTACAGTGCTAAAAAGTTATAATCCAATTGTGTTCTCTGGTGACATTAACTGAAATACATGTTCTTTTCTAGAGACCTTATTCTGTAAATCATAATAAAAAATTTCTTATGAAATTAAGTCTTTATAGATACCTCATTTAATTCTTACCATTCTATTAATAAAATACATTTTTATATTACTTATATGTATTGAGTATCTTGTCTGAAATCACACAACTAGAAAGTAGAAAGTTGTCATTTGAACACAACTCTGTGTGATTCTTAACCTTGTTTTCTTAGCCATTTTGCTGAATGGCCTTTCCTTCTTTCCTCTATCCCCCATCCCCCATCCCCCATTTCCATGTATGTGTCTGTTGTCATTCTCTACTAGCTTTTTCTTCTTATGTACAAATAGGTCAGATGTCTACATTTAAATAAACAGCATCTTTAACCCTGCTGCTCCCATACTCAAGACAAATTGACAAGTATGTTATTATCTGTGCTGTGAAGTGAAATTAGATCTTTTGATTTCAACATCAAAACCTTATAAGATGTTTTATCTTTTTTAGAAGTCCTTAATCACAAGTTAGAGCAATTAACACCTATATTACTAATAATATTTCCAGTAGTATGGATGCTTTGGTTATGGTTAAATGGTTCTGAAATCTACACCTGAGATGCTTCCCTTAAAGGGGCTCTTTGGCCAGTTCACCAATTACTTGGAAGAGGAAGAGGAAGTCTAGGTTGATTATATTGCTTTGAGTTAGAATATCTGAGAAAGTGGTACTTATATTAGTTTGCTAGTGCTGCGATGATAAAGTACTACATACTGGATGGCTTCAATAACAGAACTTTATTTTCTCACAGTCCTAGAGGCTGGAAGATCGTCTGGGACAAGTCCAGGTGCTGGCAGGGTTGATTTCTCCTGAGGCCTCTCTCTGTGGCTTGCAAGGTGGCTGTGTCTTGACATGGCCTCTTCTCTGTGTGTCCCTCCCCCTTCTTATAAGGATGTCAGTCAGATTGGATTAGGGCCCACCCTAACGGCAGCATTTTAACATTCACTTATTTAAAGACCTTATCTCTAAGAACATTACATTCTGAGGACCTGAGACTTTCAGTGTATAAATTTGGGAGGACACAATTTAGCCTATAATAGTGCTATAAAGAGAAGTTCAAAAATGAAGATGGAGATTTTATTTTTTAAG";
    std::string query = "TATTATAGTCTTCATTCTGTGTATTAGATTACTAAAGCATATTACTTCTGTCTAAATGAAATTTTGTATCCTTTAGTCAGCATCTTCCCTTTCTCCATCCACTTTTCTCTCCCAGCCTCTGGTAATCAACATTCTACCACAATTCTGTGAGTTCTACTTTTTTAGATTCCCCATGTAAGTAAGATCATGCTGTATTTGTCTTTGTGTGCCTGGCTTATCTGAGTTAACGTAATGTCCTCCAGGTTCATTCATGTTGCAAATGATAGAATTTCCTTCTTTATCAAGACTGGATAGTATTCCATTGGTGTATATATACTATATTTACTTTATGCATTCGTCTAGACACCTAGATTGCTTCCAAATCTTAGCTATTGTGAATAGCACTGCAGTTAACATGGGAGTGCAGATATTTTTTTAGCATACTGATTTCAATCCTTTGGCCCAGAAGTGGGATTACTAGATTATATGGTAGTTCTATTTTTAGTTTTTTGAGAGACCTTCATACTGTTCTCCATAATAGCTGTGATAATTTACATTCCCACCAACAATGTACAAGTTCCCTTTTTGCTACACACTCACCAACGCTTGTTATCTTTCATCTTTTTGATAGTCTTTCTAACAGGTGTGAGGTGATCCCTCATTGTGGCTTTAATCTGCATTTCCCTGGTGAGTGGTGATGATGAGCATTTTAATATATATCTGTTGGCCATTTGTACATCTTTTGAGCAATGTTTAGGTCCTCTGCACATTTTTAAATTGGGTTCTTTGTTTTCTTGCTATTGAGTTGAGTTCTTTGTCTATTTTAGATATTAGCGCCTTATCAGATATGTAGTTTGTAGATATTTTCTCTCAATCCATGGCACATCTTTTGCTCTGTTGTTTCTTTGCTGTGTATAATCTTTTCAGTTAGATGCAATCTCATATGTTTTTGCTTTTGTTGTCAGTGCTTTTAGGGTAATATTTAAGAAATCTTTGCCCAGACCAGTGTTATGGAGCATTTTCCCTATGTTTTATTTCAGTAGCTTTACAGTTTCAGGTTTTATGTTTAAGTCATTAGTCCATTTTTAGTTGATTTTGGTGTATGGTGTGAGATAAGGGTCTAATTTCTTTCTTTTGCATGTGGTTAACCTGTTTTCCCAGCACAATTTATTGAGGATTGTCCTCTTTCCATTGTATATTCTTGGCACGTTTGTTGTAAATAATTGACCACAAATGTGTGGGTTTACTTCTGGGCTCTCTACCCTGTTTGATTGGTTAGTTGGTCGGTTTTTATGCTGTGCTGTTTTGGTTACATTAGCTTTGTAACAGGTTTTAAAATCAAGTACTGTGATACCTCTAGGTTTGTTCTTTTCGCTTTGGCCATTTGGGGTTTTTTGTGGTTCCATATGAGTTTTAGGATTGTTTTATTCTGTGAAGAATGACATTGGAATTTCGTTAGGCATTGCATTGAATCTGAATATCACTTTGGGAAGTATGAATAGTTTAACAATATTCTTGCATTCCATGAATATGGATAGCTTTCCATGTGTTTGTGTCATCTACTATATCTTTCATCAGTGTTTTGTATTTTCTAGTATATACATCTATTGCCTCCTTAAATTTACTTCTAAGTTGCTTTTTCTCGATGCTACTGTAAACGGGATTGAGTTCTTAATTTCTTTTTCAGGTACTTCATTCTTAGAATAGAGAAACACTGTAGATTTTGTATCCTGCAACTTTACTGAATTTGTTTATCAGACTTAATGGTTTTTTGGTGAAGTTGTTAGCATTTTCTATATGTGAGACCATGACATCAGCAAACAGATGATTTCTCTTCTTCCTCTGATATGGATGCCTTTATTTCTTTCTCTTGCCTAATCGCTCTTGCCAGGACATCTAACTCTGTTGAATAGCAGTGGCAAGAGCGGGCATTCTTATCTTGTTTTTGATCATAGAGGAAAAGCTTTTACCTTTTTCCAGGGAGCATGGCAGCTGTGGGCTTGTTACATATGGCCTTTATTGTGTATAAACACAATATATTCCTTCTATACCTAATGTGTTGACAGTTTTTATCATGAAATAATTTTGAGTTGTGTCACATGCTTTTTCCATATCTATTTTTTTTTATTATTATACTTTAAGTTTTAGGGTACATGTGCAAAACGTGCAGGTTTGTTACATATGTATACATGTGCCATGTCGGTGTGCTGCACCCATTAACTCGTCATTTAACATCAGGTATATCTCCTAATGCTATCCCTCCCCACCTCCCCCACCCCACGACAGGCGCCGGTGTGTGATGTTCCCCTTCCTGTGTCCATCTGTTCTCATTGTTCACTTCCCACCTGTGAGTGAGAACATGCGCTGTTTGGTTTTTTGTCCTTGTGATAGTTTGCTGAGAATGATGGTTTCCAGCTTCATCCATGTCCCTACAAAGGACATGAACTCATCATTTTTTATGGCTGCATAGTATTCCATGGTGAATATGTGCCACATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGCTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAATAATCGTACGTGGGCATGTGTCTTCATAGCAGCATATTTTGTAATCCTTTTGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACACTCCCACCAACAGTGTAAAAGTGTTCCTGTTTCTCCACATCCTCTCCAGTACCTGTTGTTTCCTGACTTTTTAATGATTGCCATTCTAACTGGTGTGAGATGGTATCTCACCGTGGTTTTGATTTGCATTTCTCTGATGGCCAGTGATGATGAGCATCTTTTCACGTGTCTTTTGGCGGTATAAATGTCTTCTTTTGAGAAGTGTCTGTTCATATCCTTTGCCCACTTTTTGATGGGGTTGTTTTTTTCTTGTAAATTTGTTTGAGTTCATTGTAGATTCTGGATATTAGCCCTTTGTCAGATGAGTAGATTGCAAAAATTTTCTCCCATTCTTTAGGTTGCTGTTCACTCTGATGGTAGTTTCTTTTGTTGTGCAGAAGCTCTTTAGTTTAATTAGATGGCGTTTGTCAATTTTGGCTTTTGTTGCCATTGCTTTTGGTGTTTTAGACATGAAGTCCTTGCCCATGCCTGTGTCCTGAATGGTATTGCCTAGGTTTTCTTCTAGGGTTTTTATGGTTTTAGGTCTAACATGTAAGTCTTTAATCCATCTTGAATTAATTTTTGTATAAGGTGTAAGGAAGGGATCTGGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACCATTTGTAAAATAGGGAATCCTTTCCCCATTTCTTGTTTTTGTCAGGTTTGTCAAAGATCAGATAGTTGTAGATATGCGGCATTATTTCTGAGGGCTCTGTTCTGTTCCATTGGTCTATATCTCTGTTTTGGTACCAGTACCATGGTGTTTTGGTTACTGTAGCCTTGTAGTATAGTTTGAAGTCAGGTAGCGTGATGCCTCCAGCTTTGTTCTTATGGCTTAGGATTGACTTGGCAATGCAGGCTCCTTTTTGGTTCCATATGAACTTTAAAGTAGTTTTTTCCAATTCTGTGCAGAAAGTCATTGGTAGCTTGATGGGGATGGCATTGAATCTATAAATTACCTTGGGCAGTGTGGCCATTTTCACGATATTGACTCTTCCTACCCATGAGCATGGAATGTTCTTCCATTTGTTTGTATCCTCTTTTATTTCATTGAGCAGTGGTTTGTAGTTCTCCTTGAAGAGGTCCTTCACGTCTCTTGTAAGTTGGAATCCTAGGTATTTTATTCTCTTTGAAGCAATTGTGAATGGGAGTTCACTCATGATTTGGCTCTCTGTCTGTTATTGGTGTATAAGAATGCTTGTGATTTTTGCACATTGATTTTGTATCCTGAGACTTTGCTGAAGTTGCCTGTCAGCTTAAGGAGATTTTGGGCTGAGGCAATGGGGTTTTCTAGATATACAATCATGTCATCTGCAAACAGGGACAATTTGACTTCCTCTTTTCCTAATTGAATACCCTTTATTTCCTTCTCCTGCCTGATTGCCCTGGCCAGAACTTCCAACACTATGTTGAATAGGAGTGGTGAGAGAGGGCATCCCTGTCTTGTGCCCGTTTTCAAAGGGAATGCTTCCAGTTTTTGCCCATTCAGTATGATATTGGCTGTGGGTTTGTCATAAATAGCTCTTATTATTTTGAGATACGTCCCATCAATACCTAATTTATTGAGAGTTTTTAGCATGAAGGGTTGTTGAATTTTGTCAAAGGCCTTTTCTGCATCTATTGAGATAATTATGTGGTTTTACTCGTTGGTTCTGTTTATATGCTGGATTACATTTATTGATTTCTGTATGTTGAACCAGCCTTGCATCCCAGGGATGAAGCCCACTTGATCATGGTGGATAAGCTTTTTGATGTGCTGCTGGATTCGGTGTGCCAGTGTTTTATTGAGGATTTTTGCATCGATGTTCATCAGGGATATTGGTCTAAAATTCTCTTTTTTGGTTGTGTCTCTGCCAGGCTTTGGTATCAGGATGATGCTGGCCTCATAAAATGCGTTAGGGAGGATTTCCTCTTTTTCTATTCATTGGAATAGTTTCAGAAGGAATGGTACCAGCTCCTCTTTGTACCTCTGGTAGAATTTGGCTGTGAATCCGTCTGGTCCTGGACTTTTTTTGTTTGGTAAGCTATTAATTATTGCCTCAATTTCAGAGCCTGTTATTGGTCCATTCAGAGATTCAGCTTCTTGCTGGTTTAGTCTTGGGAGGGTGTATGTGTCCAGGAATTTACCCGTTTCTCCTAGATTTTCTAATTTATTTGTGTAAAGGTGTTTATAGTATTCTCTGATGGTAGTTTGTATTTCTGTGGGATTGGTGGTGATATCCCCTTTATCATTTTTTATTCCGTCTATCTGATTCTTCTCTCTTTTCTTCTTTATTAGTCTTGCTAGCAGTCTATCAGTTTTGTTGATCTTTTCAAAAAGCCAGCTCCTGGATTCATTGATTTTTTTGAAGGGTTTTTTGTGTCTCTATTTCCTTCAGTTCTACTCTGATCTTATTTCTTGCCTCTGCTAGCTTTTGAATGTGTTTGCTCTTGCTTCTCCAGTTCTTTTAATTGGGATGTTAGGGTGTCAATTTTAGATCTTTCCTGCTTTCTCTTGTGGGCATTTAGTGCTACAAATTTCCCTCTACACACTGCTTTGAATGTGTCCCAGAGATTCTGGTATGTTGTGTCTTTGTTCTTGTTGGTTTCAAAGAACATCTTTATTTCTGCCTTCATTTCGTTATGTACCCAGCAGTCATTCAGGAGCAGATTGTCCAGTTTCCATGTAGTTGAGCGGTTTTGAGTGAGTTTCTTAATCCTGAGTTCTAGTTTGATTGCACTGTGGTCTGAGAGACAGTTTGTTATAATTTCTGTTCTTTTACATATGCTGAGGAGTGCTTTACTTCCAAATATGTGGTCAATTTGGAATAGGTGTGGTGTGGTACTGAGAAGAATGTGTATTCTGTTGATTGGGGGTGGAGAGTTCTGTAGATGTCTATTAGGTCCGCTTGGTCCAGAGCTGAGTTCAGTGCCTGGATATCCTTGTTAACTTTCTGTCTCACTGATCTGTCTAATGTTGATAGTGGGGTGTTAAAGTCTCCCATTATTATTGTGTGGGAGTCTAAGTATCTTTGTAGGTCTCTAAGGACTTGCTTTATGAATCTGGGTGCTCCCGTATTGGGTGCATATATATTTAGGATAGTTAGCTCTTCTTGTTGAATGATCCCTTTACCATTATGAAATGGCCTTCTTTGTCTCTTCTGATCTTTGTTGGCTTAAAGTCTGTTTTATCAGAGACTAGGATTGAAACCCCTGCCTTTTTTTGTTTTCCATTTCCTTGTAGATCTTCCTTCATCCTTTATTTTGAGCCTATGTGTGTCTCTGCACGTGAGATGGGTTTCCTGAATACCGCACACTGATGGGTCTTGACTCTTTATCCAATTTGCCAGTCTGTGTCTTTTTATTGGAGCATTTAGCCCATTTACATTTAAGGTTACTACTGTTATGTGTGAATTTGATCCTGTCATTATGATGTTAGCTGGTTATTTTGCTCATTAGTTGATGCAGTTTCTTCCTAGCCTTGATGGTCTTTGCAGTTTGGCATGTTTTTGCAGTGGCTGGTACTGGTTGTTCCTTTCCATGTTCAGTGCTTCCTTCAGGAGCTCTTTTAGGGCAGGCCTGGTGGTGACCAAATCTCTCAGCATTTGCTTATCTGTAAAGGATTTTATTTCTCCTTCACTCATGAAGCTTAGTTTGGCTGGATATGAAATTCTGGGTTGAAAATTCTTGTCTTTAAGAATGTTGAATATTGGCCCCCACTCTCTTCTGGCTTGTAGAGTTTCTGCCGAGAGATCCGCTGTTAGTCTGATGGGCTTCCCTTTGTGGGTAACCCGACCTTTCTCTCTGGCTGCCCTTAACATTTTTTCCTTCATTTCAACTTTGGTGAATCTGAGAATTACGTGTCTTGGAGTTGCTCTTCTCGACGAGTATCTTTGTGGCGTTCTCTGTATTTCCTGAATTTGAATGTTGGCCTGCCTTGCTAGATTGGGGAAGTTCTCCTGGATAATATCCTGCAGAGTGTTTTCCAACTTGGTTCCATTCTCCCTGTCACTTTCAGGTACACCAATGAGACGTAGATTTGGTCTTTTCACATAGTCCTATATTTCTTGGAGGCTTTGTTCGTTTCTTTTAATTCTTTTTTCTCTAAACTTCTCTTTTCACTTCAGTTCATTCATTTCATCTTCCATCACTGATACCCTTTCTTCCAGTTGATCGAATCGGCTACTGAGGCTTGTGCATTCGTCACGTAGTTCTCGTGCTGTGGTTTTCAGCTCCATGAGGTCCTTTAAGGACTTCTCTGCATTGGTTATTCTAGTTAGCCATTCGTCTAATTTTTTTCAAGATTTTTAACTTCTTTGCCATGGGTTCGAAGTTCCTCCTTTAGCTCGGAGTAGTTTGATCGTCTGAAGCCTTCTTCTCTCAACTCATCAAAGTCATTCTCCGTCCAGCTTTGTTCCGTTGTTGGTGAGGAGCTGTGTTCCTTTGGAGGAGGAGAGGCACTCTGATTTTTAGAGTTTCCAGTTTTTCTGCTCTGTTTTTTCCCCATCTTTGTAGTTTTATGTACCTTTGGTCTTTGATGATGGTGACGTACAGATGGGGTTTTGGTGTGCATGTCCTTTGTGTTTGTTACTTTTCCTTCTAACAGTCAGGACCCTCAGCTGCAGGTCTGTTGGAGTTTGGTGGAGGTCCACTCCAGACCTTGTTTGCCTGGGTATCAGCAGCGGAGGCTGCAGAACAGCCAATATTGGTTTATACCAAATGTTGCTGCCTGATCATTCCTCTGGAAGTTTTTGTCTCAGAGGAGTACCTGGCCATGTGAGGTGTCAGTCTGCCCCTACTGGGGGGTCCCTCCCAGTTAGGCTACTCGAGGGTCAGGGACCCACTTGAGGAGGCAGTGTGTCCGTTCTCAGATCTCCAGCTGCATGCTGGGAGAACCACTACTCTCTTCAAAGCTGTCAGACAGGTACATTTAATTCTGCAGAGATTTCTGCTGCGTTTTGTTTGGTTATGCCCTGCCCCCAGAGGTGAGTCTACAGAGGCAGGCAGGCCTCCTTGAGCTGTGGTGCGCTCCACCCAGTTCCAGCTTCCCTGTGGCTTTGTTTACCTACTCAAGCCTTGGCAATGGTGGCTGCCCCTCCCGCAGCCTGGCTGCTGCCTTGCAGTTTGATCTCAGACTGCTGTGCTAGCAATGAGCCAGGCTCCGTGGGCGTAGGACTCTCCGAACCAGTTGTGGGATATAATCTCCTGGTGTGCCATTTGCTAAGACCGTTAGAAAAGCGCAGTATTAGGGTTGGAGTGACCCGATTTTCCAGGTGCTGTCTGTCACTCCTTTCTGTGACTAGGAAAGGGAATTCCATGACCCCTTGCACTTCCCGGGTGAGGTGATGCCTCGCCCTGCTTCGGCTTACGCTGGGTGCACTGCACCCAGTGTGCTGCACCCTCTGTCCAACACTCCCCAGTGAGATGAACCCAATACCTCAGTTGGAAATGCAGAAATCACCCGTCTTCTGCGTCACTCACGCTGGGAGCTGTAGACTGGAGCTGTTCCTGTTTGGCCATCTTGGCTCCACCCCTTTCCATATCTATTCAGGTGATCATATGGATTTTGTCTTCCATTCTATTAACGTGGTGTATCACATTTATTTATTTGCATATGTTGAACCATCTTTGCATCCCAAGGGTAAATCCCACTTGATCATGTGAATGATCCTTTTAATGTGCTGTTGAATTCAGTTTATTAGTATATTGGTGATTTTTGCATCCTTGTTCATCAAAGATGACATGGTAGCTTGGATTACAGGCGCCCGTCACCATGCCCAACTAATTTTTGTATTTTTAGTAGAGACAGGGTTTCACCGTATTGGCCAGGCTGGTCTTGAACTCCTGACCTCAGGTGATCCACCCACCTCGGCCTCCCAAAGTGCTGAGATTACAGGTGTGAGCCACTGCACCTGGCCTGCGGTATTTTCTTTTAGGACATTTAACAGATGTATTTGATTACTGATAAGACTTTCCATATTAAGTTGATAGTACTTTTTATTCTTGTCTGGGTTTTTATAATTTATTTTATTGGTATCTTCCAATAAGCATTTGCTATGGTTTGAATGTCCTCTCCAAAACTCTTGTTGAAATTTAATTACCATGGTGACAGTATTAAAATAGGGGACTAATATAGTTTGAGTCTGTGTCCCCGCCCAAATCTCATGTTGAATTGTAATTTCCATTGTTGGAAATGGGGCCTGGTGGAAGGTGACTGGATTATGGGAGTGGATTTCTCATGAATGGTTTAGCACCAAAACCCTTGGTACTGTCCTCACAATAGTGAGTGAGTTCTTATGTGAGCTGATAATTTAAAAGTGTGTAGTGGTCAGGCACTGTGTCTCCTGCCTGTAACCCCAGCACTTTGGGAGGCCAAGATGGGAGGATCAGTTGAGCCCAGGAGTTTGAGACCGACCTGGCCAACATAGCAAGACCTCGTCTCTATTAAAAAAAATTAAAATAAATAAATACGTAAAAGTGTGTAGTACTGCCCACTTTGCTCTCTTGCTTTCACCATGTGAAGTGTCTGCTTCCTCCTCCCCCTGCTTTCACCATGTGAAGTGCCTGCACCCACTTTGCCTCCCACCATGAGTGAAAGCTCCCTGAGGTGTCCCCAGAAGCATATGCCACCATTCTTCCTGTACAGCCTGCAGAACTGTGAGCCAGTTAAACCTCGTCTTTTATAAATTACCCAGTCTCAGGTATATCTTTATACCAATGTGAGAATGGTATACATTAAAAGAATGGTATATATTAAAAGAACCTTTTAATGGTTCCTGTATTAGGATTAGGCCATGAGGGCCTTGCCCTCATGAATAGATTAATATCATTATCACAGTAGTGGGTTTGTTACCACAAGAATGTGTTGTTATAAAAATTTTCTGTCCCTTGCTCTCACCCTCACTTGTCGTTCTGCTTTCTGCCATGGGAAGGTGTAGGACAAAGGTCCTCACCAGGAGCCAGCACCTTGATATTGGAATTTTCAGCCTCCATAACTATGAGAAATAAACTTCTTTTCTTTGTAAGTTACCTATCTTGTGGTATTCTGTTAAAGCAAAGCAAAGTGGACTAAGACACCATACCTACAGAGGGCAGCAGAGACATAGCAGTGAATTCTCATAACCTAATTTTTCAGGAAGTGTGAAGGAAATTGGATTTAGGTGAACTTTGTTATGTCATTCTTATATACTTTATGTAGTGACTCCTATGCAAATACAGGCTAATAGACATAATCAATGTAAAGTAATGTGAGTGAGACTCAGGCACATTTGTGAGAAGGTAAACAACAAGCTGGCTGGCACAGATAAGAGATCTGAAAGTAGTTCATGTAGTTTTAAAGGAAAAGTAAATGTGGCTTTCAAATTTGAAATACGGATTCATTTTTATGAAACAGATATTTTGTTTGTATTCTTATCAGTAAGGGTAGAAAAAATGATTCAGTAACTCATAGAAGATGGAATTACGGAAAAGTTTAAAAAAACAGGAATAAACAAAATTTTGTACTACATGGGAGTGTATTATAATCTTGCAGAAAGTGATCTGACAAACCTTTAAAAACATTTCTAAATTTAATGTGCATATGTATAATTTTCTGTATTTTTACTTATTTATATAAATTTTTGAAAATCAGTTTTGCTTTTAATCGGTTTTAATAGTGAGTTAGAATAAAATGATATGACTAGATGTGTACTAATTGTTAATGCCTAAGAACTGCCCATTGTTTTCATCTTTTTTCCCAATGTTTTACATTTTAAAATACACATAACATAAAATTTACTATTTTGACCATTTTAAAATGTGCAATTCAGTGGTATTGAGTATCTTCATAATGTTGTGTAACTACCATCACCACCTTCTATCTCCAGAACTCTTTTCATCTTGTAAAATTGAAATTCTGTACTCATTAAACAGTAACTCCTCCTTCTGTCTTTCCCTAGCTCCTGGCAACCACCATTCTACTTTGTGTATGATTTTGACTACTGTAATTACCTCATATGAATAGAATCATACAGTATTTGTTTTTCTGTGACTGGTTTATTTCACTTAGCATAATATCCTTAAGGTTAATCCATGTTGTTAGAATTTTTTATTTTTATTTTTTGAGATGGAGTCTCACTCTGTTGGCCAGGCTGGAGTGCAGTGGTGCGATCTTGACTCACTGCAACCTCCGCCTCCTGGGTTCAAGTGATTGTTGTGCCTTGGCCTCCTAGTATCTGGGATTATAGGCGCACGCCACCACGCCTAGCTAATTTTTGTATTTTTACTAGAGACAGGGCTTCACTATGTTGGCCAGGCTGGTCTTGAACTCTTGACCTCAAGTGATCCGCCCACCTCAGCCTCACAAAGTGCTGGGGTTACAGGCATGAGCCACCGAGCCCGGCTAGAATTTCCTTCTTTTTAAAGGCTGAATAATATTCCATTGTATGTTTATACCACATTTTTTACTTAATTTGTGACCTAACATATGGTCTATCCTGGAAAATGTTCTATATTCACTTGAGAAGAATATGTATGCTATTGTTGGGTACAGTGTTGTGTATATGTCTGTCAGATCTAGTTGGTTATTGTCTGAGGTCCTCTGTTTCCTTACTTATCTTCGATCTGGTTATTCTATCCACTATTGAAAGTAGGGTATTGATTTAATCATCTGTTATTTAGAAAGGTCTGTTTCTCCCTTCAGTTGTCTGTTTTTCCTTCATGTATTTTTGAGGTCTTTTAATTGGTGCATAAATGTTTATAATTTTTGTCTTCTTGTTGTATCAAATCATGTATTCATATATAATGTCCTTTGTCTCTTGTAACCTTTTTTGATTTAGACTGTTTTCTCTGTTGTTAGTATAGCCACCTCTACTCTGTTCTGGTTAGTATTTGCATGGAATATCTTTTTCTGTCTTTTAATTTTCAAACTGTGTCTGTCAGTCTAAAGTAAGTCTAGTATAGGCAACATATGGTTGAATCCTCCTCCATTCTGTCAATCTCTGTCTTTTGATTGGAGAGTTTAATCCATTTGAATTTAAATTATTTGTACAGAAGGATTTCTGTCATTTTGTTGTTTGTTGTTTACATAACTTACTGCTTTTTGTGGGGGCCCTCATTTCCTGCATTACTGTCTTCTTTTGTGTTTACTTGATTTTTTTCATGAAATGTTTAAATTCCTTTCTTCTTTCTTTCTTCTTTTTTTTTTTTAATGGAGTCTCTATCGCCCAGGCTGGAGTGCAGTGGCATGATCTTGGCTCACTGCAACCTCTGCCTCCTGGGTTCAAGCAGTTATCCTGCCTCAGCCTCCTGAGGCTGGGATTACAGGCATGTGCCACCACGCCTGGCTCATTTTTGTATGTTTAGTAGAGATGGGGTTCCACAATGTTGGCCCAGCTGGACTCGAACTCCTGACCTCAGGTGATCCACCCACCTTGGCCTCCCAAAATGCTGGGATTACAGGTGTGAGCCACCGCACCCAGCCTCTTTTCTCATTTCGTTTTTTTTTTTTTTTTTTTTTTTGAGACAGAGTCTCACTCTGTCACCCGAGCTGGAGTGCAGTGGCACAGTCTCAACTCACTACAACCTCTGCCTCCCGGGTTCAAGCGATTCTCATGCCTCAGCCTCCCTAGTAGCTGAGATTACAGGCGTGTGCCACCATGCCCGGCTAATTTTTGTATTTTTAGTAGAGATGGGGTTTTGCCATGTTGGCCAGGCTGGTCTCAAACTGCCAACCTCAGGTGATCTTCCCGCCTTGGCCTTCCATAGTGCTGGGATTACAGTCATGAGCCACTGTGCCCGCCCTAGTACTTTGAACAACTTTAAAACAGTTGTCATAAAGTCTTTCTCTGGAAGGTCTGCCATCAGGTCTTTTTCAAGGACAGTTGTGTTGTTTTATTTTTTTTTCTTTTGAATGTGTCATACTTTCCTGTATATCTTGTTATTTTTTTGTTGAAAACAGGGCATTTGAATATAATGATGTGGTAAGTCTGGAAATTAGATTCTCCCTCTTCCCCAGCATTTGCTGTTTGTTTTGTTTTGTTTTGAGGCAGAGTCTCACTCTGTTGCCCAGGCTGGAGTGCAGTGGTGTGATCTCGGCTCACTGCAACCTCTACCTCCTGGGTTCAGCCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGATTACAGGTGTGTGCCACCATGCCCAGCTAATTTTTGTGCTTCTTAGTAGAGATGGGGTTTCGCCATGTTGCCCAGGCTCGTCTCAAACTCCTGACCTCAGGTGATCTGCCCACCTCAGCCTCCCAAAGTGCTGAGATTATAGGCATGAGCTACCATGCCCGGCCTGCTGTTTTTTTTTTTTTTTTCTAATGTTACTGTTTTTGTTTATTTTTATTGTTTGAGGCTTGTCTGTATGTGAAGGATCAGCCTGAGGTGTAAGCTTAAGATTCTCTGAGTCTTTTCTAAGCCTGAGTTTTTCTCCAGGAATGTGCAGGAAGTTTCTAATTTTTCCCTTATATACAGTTGCTTTCAAATGCCGTAGTCTTTAACATCTGACTCCCAAAAGGGAAAAAGAGAAAAACTAAATCCTTTAAATCCCCAGCAGTTACTTTGTGGGTGGGGGTGGGGATGGGGGTGGGAGAGGGCTTGTAACAGTGAGTGGAGATGAAGCAACAATGTGTGTCCACCTCTTTGTCCCATTTCTGTGATCAGAAGTAGCATTCAGTGATCAGGGCACACATGCCCTATATTTGGAAGAGAGGGTCCTTTTTGCCTGTCCTTGGTTCCTGTAAGCTGCGCGCAAGCTGCTCCTTGGAATGCATGCACAACTGCCTGCCTTGGGGCTGAGGGTGGCAAATGGATAGCTGCTACTGTGCTAGTGCTAATACCTGCAAGTGATTAAAATTATATTCAGGCCTTCCCCTGGACATTAACAAGCCTTCAGTAGATTCCAGAGTTTCAAAATAGTTAAATAGAACAGATTCTGCCAGTGCAATGATCTTCTAGGTGGAGCAACACATGCGTGGTGCTTCCTACTCCACCATCTTCCCAGAATCCTCTGATTGCGTATAAGTTACCTATAGTTGAAAGAGTTTAATTCTGATCATCATACATTTTAGAAACCTCTCATATATGTTTATATAAAAACTGATAAAATTTAGTGCTTGTAAAGGGAAATCTGGAGTTCTTAGTCAAGCTTTTTGTATGATTATATTTTAAAGATTTAATTATAGAAACATTGTACTGAAGTCTATATTCCAGTGTGTTCTTATGGTAAAAAGAAAAGGTTGGAGATTTTATTGAGGAGTCTGTATTTTTTTACTGTTATTTTATATTATAACTCATGTATTATATACTGTATCTGTTCATTATCTTAGTATAGCATTTTCGGGCTTTTATTCTCTTTATGAAAGTGAGTAACTGTAGGTACCAATATCAATTGTTGCATCTGTGCCAGAAATTATATTTTATCTTTTCATCATCTCCCTTCTTGTATAGCTGAACTTTAGCTTAAAAGATGCCAGGACAAGGAGTTTGTCAGATTGTTCCTTTTCTTTGCATGAACCTTCTTCCCCTTTCAGTATCACTGGGTCATCTTATTTCCCTTCATTCCCCTCTTTCTTCTCCATCTTTTCATTTCCCAAATACTTCTGTTTTCTGATATAAGAAGGACTATAACAGAAAATGATTGCTTTAGGAAATATTGTTACAACTGTGTAATGTTAAAAAGTAGTAAGTTGGCCAGGCACAGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCTGAAGTGGGCGGATCACGAGGTCAGGAGGTCAAGACCATTCTGGCTAACACAGTGAAACCCTGTCTCTAAAGAATACCAAAAAAAATTAGCCGGACATGGTGGCGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCAGAGCTTGCTGTGAGCCGGGATTGCACCACTGCACTCCAGCCTGGGCGACAGACTCTGTCTCAAAAAACAAAGTAATAAGTCGGCCAGCGCTGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCAGGCAGATCACGAGGTCAGGAGATTTGAGACCATCCTGGCTAACACAGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGTTGGTGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGAAGAATGGTGTGAACCCAGGAGGTGGAGGTCGCAGTGAGCTGAGATCGCGCCGCTGCACTCCAGCCTGGGTGACAGAACGAGACTCCATCTCAAAAAAAAAAAAAAAAAAAAGTCACTATGGATGTCACTTTTGTAGTCTTCTGGTATTATAAACAATTCTTTGGAGTTGTGCTTTTAATAATTACAGAGGGGCACATGTCGGACTACCTTTCCTGAAAATAACAATGACAACTTTTGGACAAAATATAAAAAAAGTATGTGAAAGCACAGTAGGAGGACAGAAAGTAGGAAGAAACTGAAGGATGTTGACATTTGGAAGATAGAATGGTATAAGATGAGCATTTCATTTTTTGTGGATGTTTCCCAGAGGACAGAACTCGGTCATTGGTTGTGTAGGTGACGGAAGCTTGAATAAAAACCACAAGACAGAGTTAGGGGTGACTCTAGCAGCTAGAAAGTTAGTGAGTATAAAGTAGTCTTTCATAAAGAAGGGTGTCCCAAATTCTGCATAAAAGTCTTAGGCTAAAATAACTAAAGATATTTCATTTGCACAAGAAAGAGTTCAGAGGTGTATAGCCAACTTACCTACCTATTAAGCAAAAAAAAACTAAATAGTCTTCAGAAGAACATAATCAAATCCAGAGTTTCTAAAAATGTACTATGCACAGTGTCCAGAATACAATCCATATTCAGTAGACACACACATACACAAACAGAAAAATGTGACCCATACTTAAGTGAAAAGCCAATTAATTGAGACTGGCCTTGATATGACTCCTATGTTAGAATTAGCAGCTGTTATAACTCTTCTTAAGGTTATATAAAGAAAATATGTACACCATGAATAAACAGATAAGAAACCTCAGCTGAAAAATGGAAACTAGGTAAAAGGACTAAATGAAAATTCTATAGGTAAAAAATACAATATTTCAAATGAAAACTGATTAGAGATTACAGAAGATTAGGAAAAGAAAAGGAGTCATTGAAGTTGTTGGATTCATAGAAGTTATTCTTTCTTAAAAAGTGAAAAAGAGATTTCAAAAAATGAAGCTAGCCTCAGTGGCTGGTAGGACAATATCAAAAGGTCTAACACTTGGAGTCCTAGAAATAGAGGCCAGGGGAAATTGGGTCAAAAAATTCTGGAGAAATTATGGCTGGAAATTTCCCTAAATTGGTAAAAAACATGAATTTACATGTTGAAGAAGTTCAGTGAATTCCAAACAAGATAAATACAATAAAATAACACTTAAGTACATTATAGTGAAGCTGTTGAAAATAATAGACAAAAAATCCTGCCTAGTTGCAGTGGCTTATGCCTGTAATGCCAACACTTGGGGAGGCCAAAGAGGGGTATTGCTTGAGCTAAGGAGTTTGAGACCAGCCTGGGAAACGAAGTGAGACTTCATCTCTACAAAAAGAAAAAAATCTTGAAAGCACCCTGAGAAATGCTAATGAATGGTAGCTGACTTCTTTTTTTTTTTATATTATTATTATACTTTAAGTTTTAGGGTACATGTGCACAATGTGCGGGTTTGTTACATATGTATACATGTGCCATGTTGGTGTGCTGCACCCATTAACTCATCATTTAACATTAGGTATACCTCCTAATGCTATCCCTCCCCACTTCCCCCATCCCACAACAGTCCCCGGAGTGTGATGTTCCCCTTCCTGTGTCCAAGTGTTCTCATTGTTCAGTTCCCACCTATGAGTGAGAACATGCGGTGTTTGGTTTTTTGTCCTTGCGATAGTTTGCTGAGAATGATGGTTTCCAGTTTGATCCATGTCCCTACAAAGGACATGAACTCTTCATTTTTTATGGCTGCATAGTATTCCATGGTGGATATGTGCCACATTTTCTTAATCCAGTCTATCATTGTTGGACATTTGGGTTGGTTCCAAGTCTTTGCTATTGTGAATAGTGCCGCAGTAAACATAGGTGTGCGTGTGTCTTTATAGCAGCATGATTTATAGTCCTTTGGGTATATACCCAGTAATGGGATGGCTGGGTCAAATGGTATTTCTAGTTCTAGATCCCTGAGGAATCGCCACACTGACTTCCACAATGGTTGAACTAGTTTACACTCCCACCAACAGTGTAAAAGTGTTCCTGTTTCTCCACATCCTCTCCAGCACCTGTTGTTTCCTGACTTTTTAATGATCACCATTCTAACTGGTGTGAGATGGTATCTCATTATGGTTTTGATTTGCATTTCTCTGATGGCCAGTGATGGTGAGCATTTTTTCATGTGTTTTTTGGCTGCATAAATGTCTTCTTTTGAGAAGTGTCTGTTCATATCCTTCACCCACTTTTTGATGGGGTTGTTTGTTTTTTTCTTGTAAATTTGTTTGAGTTCATTGTAGATTCTGGATATTAGCCCTTTGTCAGATGAGGTAGGTTGCGAAAATTTTCTCCCATTTTGTAGGTTGCCTGTTCACTCTGATGGTAGTTTCTTGTGCTGTGCAGAAGCTCTTTAGTTTAATTAGATCCCATTTGTCAATTTTGGCTTTTGTTGCCATTGCTTTTGGTGTTGTAGACATGAAGTCCTTCCCCATGCCTATGTCCTGAATGGTATTGCCTAGGTTTTCTTCTAGGGTTTTTATGGTTTTAGGTCTAACATTTAAGTCTTTAATCCATCTTGAATTAATTTGTGTATAAGGTGTAAGGAAGGGATCCAGTTTCAGCTTTCTACATATGGCTAGCCAGTTTTCCCAGCACCATTTATTAAATAGGGAATCCTTTCCCGAAAACTGGCACAAGACAGGGATGCCCTCTCTCACCACTCCTATTCAACATAGTGTTGGAAGTTCTGGCCAGGGCAGTCAGGCAGGAGAAGGAAATAAAGGGTATTCAATTAGGAAAAGAGGAAGTCAAATTGTCTCTGTTTGCAGATGACATGATTGCATATCTAGAAAACACCATCATCTCAGCCCAAAATCTCCTTAAGCTGATAAGCAACTTCAGCAAAGTCTCAGGATACAAAATCAATTTGCAAAAATCACAAGCATTCTTCTAGACCAATAACAGACAAACAGAGAGCCAAATCATGAGTGAACTCCCATTCACAATTGCTTCAAAGAGAATAAAATACCTAGGAATCCAAACTTACAAGGGACGTGAAGGACCTCTTCAAGGAGAACTACAAACCACTGCTCAATGAAATAAAAGAGGATACAAACAAATGGAAGAACATTCCGTGCTCATGGGTAGAAAGAATCAATGTCGTGAAAATGGCCATACTGCCCAAGGTAATTTATAGATTCAATGCCATCCCCATCAAGCTACCAATGACTCTCTTCACAGAATTGGAAAAAACTACTTTAAAGTTCATATGGAACCAAAAAAGAGCCCGCATTGCCGAGTCAATCCTAAGCCAAAAGAACAAAGCCGGAGGCATCACGCTACCTGACTTCAAACTATACTACAAGGCTACAGTAACCAAAACTGCATGGTACTGGTACCAAAACAGAGATATAGATCAGTGGAACAGAACAGAGCCCTCAGAAATAATGCCACATATCTACAACCATCTGATCTTTGACAGACCTGACTTCTTATCAGAAATAATGCAGGCCAGCAGGCTATAAACACCTTTAAAGTACTGGGAAAATAAAAATCAACTCAGAATTCTGTAGCCAAATAAGAAAACATACTTTCTAGCTAAGAAAATCTTTATACTACAAGAAATGTTTAAGGAGATACTTCAAGCTGAAGGAAAGTGTTACCAGCTGAAACTCAGAACTTCAGGAATGAAGAGCACTGCAATTGGTAAAAATGTAAGTCTTCCCCACTCCATTAAGAAATATACAGTTGTGACTGTTTAAAGCAGAAATTATAATGTTACAATATGGCATTTATATATGACTGCTGTGTATAAAGGATGAGTCTAGGTTGATGGAGCTGTGTGCTTACTGTGGTCTTTATTTTATGAAGGGGTAGAATTAGCTTTAAGTACACTGCAAAGTTAAAGATATACATTATAATTCCTAGAACAAACACTGGAAAAGAAAAAAATGCAAAGAGGTGATAACTAAAAAACAAGAGAAGTCAAAATGGAATTCTGAAATGTATTTAGAGAAAATGAAAAAAGGTATAAAGGAAGGAACAGCAGAATATAAAACAAATGAGATGTGTGTCAGTCTGTTCTCATACTGCTAATAAAGACATACTCGAAACTGGGTAATTTATGAAGGAAAGAGGTTTAATTGACCCAGTTCAGCATGGCTGGGGAGGCCTCAGGAAACTTACAATCATGGAGGAAGGGGAAGCAAACATGTCCTTCACAGGGTGGCAGCAAGGAACAGAATGAGAGCAAAGTTGGGGGAAAGCCCCTTATAAAACCATCAGATCTCCTTAGAACTCACTCATTATCACGAGAATAGCATGAGGGTAACTGCTCCATGATGCAGTTACCTCCCACCAGGTCCCTTCCACGACGTGTGAGGATTATGGGAACTACAATTCAAGATGAGATTTGGGTGGAGACACAGCCAAACCATATCAGGATGTAAAGAAAACAAATAGCAAACTTACAAACCTAGGTCAGGTCTATGAATAATTACCTTAATTGAAAATAGATTAAGTGTTTTATTTGAAGGTAGAGATTCATACATTGGATAGAAAACAAGACTTCCTGTATTCTGTCTACCAGAGATGCATTTTAAAAATAATGACACAGATAGATTGAAAGTAAATGTATGAGGAAAGGTGTAGCATGCAAACAGTAAGTATAAGATGGAGTTGTTATATTAATATTAGATAAGTAGACTTCAAGACAGTACTGGAGATAAAGAAGGGTATTTCATAGTGATAAAACTCATTAGAAATACATAGTAATCATAAATTTGTATATACCTAATAAGATCTACAGTCATAGGTGGAAATTTTAACACTACTCTGTAACTTATAGAAAACCTAGACAAAAAATCAGTTAGAATATAGAAGATCTGTAAAAATACTTCAGCCACCTTGACCTAACTGACATTTGTAGAATATCATGGGCAACTAAAATAGAGTGTACATTCTATTCAAGTACACATGATATGCATACCGAGATAAATAATATATGGGTCCTAAAACATGTCTCACTACATTTCAAAAGTTTGAAATCCTAAAGAGTATACTCTGTAACTACAGTAGAATTTAAAGTAGAAACCAATAGCAATAAGATATTCAGGAAAAATCTGAATATTTGGAAGTTAAACAGCACAGTTATCCCTAGTCGATCTATGAATTTAATGTGGTTCAAATCATTTTCTTAGCAGGACTTTTCTCTTTTGGTAGAAATTGACAAGCTGATTATAACGTTTATATGAAGATAGAAAGGACCTATGTTATACAAAGCAATCCTGATGATTTACAAGGTGGAAAACTTAAATTTCCTTATTTCAAGATTTACTCTAATAATACAGTAATGAAGACTGTGATATCAGCATAAGGATATAAAAGCATATAAATGGAGCAGAACAGATAACTGATAAATTATTTTTAACAAAGACAGCTATCCAATGAGGGAAAAGAAAGGCTTTTAAATAAAAAGTGCTGGAACAACTAGACATCTGTATGGATGAAAGTTCACCTTGACTCCCTACCTCACATGGTACCCAAAAGATAAATCAAGATGAATCATAGATCTATATGTTGAAGCTAAAGTATAAAGCTTTTAGAATAAACAAAGAACAGGCAAAGGTTTCCTTGATAGGATACAGAAAGCAGCAACTGTTTTAAAGAAAACACCTAACACTCTATAAACTTCTTCAAAATAAAAAACTTCTGCTCATTAAGAGACATGAAGAAAATAATAGCAAAACCTGTATCTGTAAAGGACTTGTATTAAGACTATGTAAAGAAATTCTACAACTGAACAATGAAAATACAACTTGATTAATTGGAGAAGGATTAGTACACACACTTCATAAATGATGCTATATTGAATGGCCAGTGAACACATGAAAATGTTACTAAACATCGTCAATTATTAAGGAAATACAAAAGTAAACTGCAAGGAGGTGCTGCTTTACTTGTCCCATTAGCTTAAAGAGTAATACTACCATATATGTAGAGCAATGAGAACTCTTATAAAATTAAGAGTTTTTCTACCCTAACCCATGACTTAGCAATTCTACTCCTAGTTATTTACTTAATAGAAAAAGAAAACACATCTACAAGATTTGTTCAAAAGTGTTTATAATAGGCCGAAACTGGACACAGCTCATATACCCACCAGTTGGAGACTGGATCAACAAATAGGTATAGTCATAAATGGAATACAACTTAGCAATGATAAGGAACACTGTCTTGGTACATCAACAACATAATGAATTTCAAATACATGGTAGGTGAAAGAAACCTTACACAAAAGAGTACATTCTGTATTCCATTTACATCGAGTTCTAGAACAGGCAAAATTGATACATAATGGAAGAAAAAGCAGTGCTGTCTTCGATAGTGGGGTGGGGATTAATTGCGAAAAGTCTTAGAGGAACTTGCTGGGTGATAGGAAGTATTTTATTTCTTGTTTGGAGTTTGAGTTCAGTGCTATATGTATATGTTAAAACTTATCAAATGGTGTGTCTTTTATTGTATATAAAACCTACTGAAAAAAGTAAATACATATTAATAGAAATATACATTGAATTCTTGTTATTGATGTTCATGCTGAGATATAATGGGATGATGTGTACTGTCATCTGCAACTTTGAAATGCATCGAAAAGTAAGATGGATTGGTGGAGCCAAGATGGCTGAATGGGAACAGCTCCAGTCTAGAGCTCCCAGCTTGAGCGACGTAGAAGACGAGTGATTTCTGCATTTCCAACTAAGGTACTGGGTTCATCTCACTGGGGAGTGTCGGAAAGTTGGTGCAGGACAGTGGGTACAGTGCACCCAGCGTGAGCCGAAGCAGGGCGAGGCATCACCTCACCCAGGAAGCGCAAGGGGTCAGGGAATTCCCTTTCCTAGATGGAGAAAGGGATGACAGATGGCACCTGGAAAATTGGGTCACTCCCACCCTAATACTGCGCTTTTCCAATGGTCTTAGCAGATGGCACACCAGGAGATTGTATCCCACGCCTGGCTCGGAGGGTCCTACACCTGCGGAGCCTCGCTCATTGCTAGCACAGCAGTCTGAGATCAAACTGCAAGGCAGCAGCGAGGCTGGGGGAGGGGCGCCCACCATTGCCAAGGCTTGAGTAGGTAAACAAAGCCACAGGGAAGCTGGAACTGGATGGAGCCCACCACAGTTCAAGGAGGCCTGCCTGCCTCTGTAGACTCCACCTCTGGAGGCAGGGAATAACCAAGCAAAAGGCAGCAGAATCCTCTGTAGACTTAAATGTCTGCAGCTTCGAAGAGAGTAGTGGTTCTCCCAGCACGCAGCTGGAGGTCTGAGAACAGAAAGACTGCCTCCTCAAGTGGGTCCCTGACCCCCCCGAGTAGCTTAACTGGGAGGCACCCCCTAGTAGGGGCAGACTAACACCTCACATGGCCGGGTACTCCTCTGAGACAAAACTTCCAGAGGAACAATCAGGCAGCAACATTTCCTGCTCACCGGTATCCGCTGTTCTGCAGCCTCCACTGCTGATACTCAGGGAAACAGGGTGTGGAGTGGACCTCCAGCAAACTCCAGCAGACCTGCAACTGAGGGTCCTGACTGTTAGAAGGAAAACTAACAGAAAGGATATCCACTCCAAAACCCCATCTGTACGTCACTATCATCAAAGACCCAAAGGTAGATAAAACCACAAAGATGGAGAAAAAACAGAGCAGAAAAACTGGAAACTCTAAAAATCAGAGTGCCTCTCCTCCTCCAAAGGAACACAGCTCCTCACCAGCAACGGAACAAAGCTGGACGGAGAGTGACTTTGACGAGTTGAGAGAAGAAGGCTTCAGACGATCAAACTACTCCGAGCTAAAGGAGGAAGTTCGAACCCATGGCAAAGAAGTTAAAAACCTTGGAAAAAAACGAGACAAATGGCTAACTAGAATAACCAATGCAGAGAAGTCCTTAAAGGACCTCATGGAGCTGAAAACCACAGCATGAGAACTACGTGACGAATGCACAAGCCTCAGTAGCCGATTCGATCAACTGGAAGAAAGGGTATCAGTGATGGAAGATGAAATGAATGAAATGAAGTGAGAAGAGAAGTTTAGAGAAAAAAGAATTAAAAGAAACGAACAAAGCCTCCAAGAAATATAGGACTATGTGAAAAGACCAAATCTACGTCTCATTGGTGTACCTGAAAGTGACAGGGAGAATGGAACCAAGTTGGAAAACACTCTGCAGGATATTATCCAGGAGAACTTCCCCAATCTAGCAAGGCAGGCCAACATTCAAATTCAGGAAATACAGAGAACGCCACAAAGATACTCGTCGAGAAGAGCAACTCCAAGACACGTAATTCTCAGATTCACCAAAGTTGAAATGAAGGAAAAAATGTTAAGGGCAGCCAGAGAGAAAGGTCGGGTTACCCACAAAGGGAAGCCCATCAGACTAACAGCGGATCTCTCGGCAGAAACTCTACAAGCCAGAAGAGAGTGGGGGCCAATATTCAACATTCTTAAAGACAAGAATTTTCAACCCAGAATTTCATATCCAGCCAAACTAAGCTTCATGAGTGAAGGAGAAATAAAATCCTTTACAGATAAGCAAATGCTGAGAGATTTGGTCACCACCAGGCCTGCCCTAAAAGAGCTCCTGAAGGAAGCACTGAACATGGAAAGGAACAACCAGTACCAGCCACTGCAAAAACATGCCAAACTGCAAAGACCATCAAGGCTAGGAAGAAACTGCATCAACTAATGAGCAAAATAACCAGCTAACATCATAATGACAGGATCAAATTCACACATAACAATATTAACCTTAAATGTAAATGGGCTAAATGCTCCAATAAAAAGACACAGACTGGCAAATTGGATAAAGAGTCAAGACCCATCAGTGTGCGGTATTCAGGAAACCCATCTCATGTGCAGAAGACACACACAGGCTCAAAATAAAGGATGGAGGAAGATCTACCAAGGAAATGGAAAACAAAAAAAGGCAGGGGTTTCAATCCTAGTCTCTGATAAAACAGACTTTAAACCAACAAAGATCAAAAGAGACAAAGCAGGCCATTATGTAATGGTAAAGGGATCAATTCAACAAGAAGAGCTAGCTAGCCTAAATATATATGCACCCAATACGGGAGCACCCAGATTCATAAAGCAAGTCCTTAGAGACCTACAAAGATACTTAGACTCCCACACAATAATAATGGGAGACTTTAACACCCCACTATCAACATTAGACAGATCAGTGAGACAGAAAGTTAACAAGGATATCCAGGTATTGTACTCAACTCTGCACCAAGCGGACCTAATAGACATCTACAGAACCCTCCACCCCAAATGAACAGAATATACATTCTTCTCAGCACCACACCGCACTTATTCCAAAATTGACCACATAGTTGGAAGTAAAGCACTCCTCAGCAAATGTTAAAAGAACAGAAATTATAACAAACTGTCTCTCAGACCACAGTGCAATCAAACTAGAACTCAGGATTAAGAAACTCACTCAAAACTGCTCAACTACATGGAAACTGAACAACCTGCTCCTGAATGACTACTGGGTACATAACGAAATGAAAGCAGAAATAAAGATGTTCTTGGAAACCAACGAGAGCAAAGACACAACATAACCAGAATCTCTGGGACACATTCAAAGCAGTGTGTAGAGGGAAATTTATAGCACTAAATGCCCACAAGAGAAAGCAGGAAAGATCTAAAATTGACACCCTAACATCACAATTAAAAGAACTAGAGAAGCAAGAGCAAACACATTCAAAAGCTAGCAGAAGGCAAGAAATAACTAAGATCAGAGCAGAACTGAAGGAAATAGAGACATAAAAAACCCTTCAAAAAATCAGTGAATCCAGGAGCTGGTCTTTTGAAAAGATCAACAAACTTGATAGGCTGCTAGAAAGATTACTAAAAGAGAAAAGAGAGAAGAATCAAATAGATGCAATAAAAAATGATAAAGGAGATATCACCACTAATCCCACAGAAATACAAACTACCATCAGAGAATACTATAAACACCTCTACGCAAATAAACTAGAAAATCTAGAAGAAATGGATACATTCCTGGACACATACACCCTCCCAAGACTAAACCAGAAAGAAGTTGAATCTCTGAACGGACCAATAACAGGCTCTGAAATTGAGGCGATAATAATAGCTTACCAAACAAAAAAAGTCCAGGACTGGATGAATTTACAGCCGAATTCTACCAGAGGTACAAGGAGGAGCTGGTACCATTCCTTCTGAAACTATTCCAATCAATAGAAAAAGAGGGAATCCTCCCTAACTCATTTTATGAGGCCAGCATCAGCCTGATACCAAAGCCTGGCAGAGACACAGCAAAAAAAGAGAATTTTAGACCAATATCCCTGATGAACATGGATGCAAAAATCCTCAATAAAATACTGGCAAACCAAATCCAGCATCACATGGAAAAGCTTATCCACCATGATCAAGTGGGCTTCATCCCTGGGATGCAAGGCAGGTTCAACATATGCAAATCAATAAACATAATCCATCATATAAACAGAACCAACGAGTAAAACCACATGATTATCTCAATAGATGCAGAAAAGGCCTTTGACAAAATTTAACAACCCTTCATGCTAAAAACTCTCAATAAATTAGGTATTGATGGGACGTATCTCAAAATAGTAAGAGATATCTATGACAAACCCACAGCCAATATCATAATGGGCAAAAACTGGAAGCATTCCCTTTGAAAACTGGCACAAGACAGGGATGCCCTCTCTCACCGCTCCTATTCAACATGGTGTTGGAAGTTCTGGCCAGGGCAATCAGGCAGGAGAAGGAAATAAAGGGTATTCAATTAGGAAAAGAGGAAGTCAAATTGTCCCTGTTTGCAGATGATACGATTGTATATTTATAAAACCCCATCGTCTCAGCCCAAAATCTCCTTAAGCTGATAGGCAACTTCAGCAAAGTCTCAGGATACAAAATCAATGTGCAAAAATCACAAGCATTCTTATACACCAATAACAGAGAAAGAGAGAGCCAAATCATGAGTGAACTCCCATTCACAATTGCTTCAAAGAGAATAAAACACCTAGGAATCCAACTTATAAGGGATTTGAAGGACCTCTTCAAGGAGAACTACAAACCACTGCTCAATGAAATAAAAGAGGATACAAACAAACAGAAGAACATTCCATGCTCATGGACAGGAAGAATCAGTATCATGAAAATGGCCATACTGCCCAAGGTAATTTATAGATTCAGTGCGATCCCCATCAAGCTACCAATGACTTTCTGCACAGAATTGGAAAAAACTACTTTAAAGTTCATATGGAACCAAAAAAGAGCCGGTATTGCCAAGTCAATCCTAAGGTGAAAGAACAAAGCTGGAGGCATCGTGCTACCTGACTTCAAACTATACTACAAGGCTGCAGTAACCAAAACTGCATGGTACTGGTACCGAAACAGAGATATAGAGCAATGGAACAGAACAGAGCCCTCAGAAATAATGCCACATATCTACAACCATCTGATCTTTGACAAACCTGACAAAAACAAGAAATGGGAAAACGATTCCCTATTTAATAAATGGTGCTGGGAAAAGTGGCTAGCCATATGTAGAAAAGTGAAACTGGATCCCTTCCTTACACCTTATACTAAAATTAATTCAAGATGGATTAAAGACTTAAATGTTAGACCTGAAACCATAAAAACCCAAGAAGAAAACCTAGGCAATACCATTCAGGACATGGGCATGGGCAAGGACTTCGTGTCTAAAACACCAAAAGCAATGGCAACGAAAGCCAAAATTGACAACTAGGATCTAATTAAACGAAAGAGCTTCTGCACAGCAAAAGAAACTACCATCTGAGTGAACAGCAACCTACAGAATGGGAGAAAATTTTTACAGTCTGCTCATCTGACAAAGGGCTAATATCCAGAATCTACAATGAGCTCCAACAAATTTTCAAGAAAAAAACGCCATCAGAAAGTGGTCGAAGGATATGAGCAGACACTTCTCAAAAGAAGACATTTATGCAGCCAAAAGACACTTGAAAAGATGCTCATCATCCCTGGCCATCAGAGAAATGCAAATCGAAACCACAATGAGATACCATCTCACACCAGTTAGAATGGCAGTCATTAAGAAGTCAGGAAACAACAGGTGCTGGAGAGGATGTGGAGAAATAGGAACACTTTTACACCGTTGGTGGGACTGTAAACTAGTTCAACGATTGTGGAAGTCAGTGTGGCGATTCCTCAGGGATCTAGAACTAGAAATACCATTTGACCCAGCCATCCCATTACTGGGTATATACCCAAAGGATTATAAATCATGCTGCTATAAAGACACATGCACACGTATGTTTATTGGGGGACTATTCACAATCGCAGAGACTTGGAGCCAAGCCAAATGTCCAACAATGATAGACTGGATTAAGAAAATGTGGCACATATACACCATGGAATACTATGCAGCCATAAAAAAGGATGAGTTCATGTCCTTTGTAGGGACATGGATGAAGGTGGAAACCATCATTCTCAGCAAACTATCAGAAGGACAAAAAACCAAACAGCGCATGTTCTCACTCATAGGTGGGAATTGAAAAGTGAGAACACAGGGGCACAAGAAGGGGAACATCATACGCCGGGGCCTGTTGTCGGGTGAGGGCAGGGTGGAGGGATAGCATTAGGAGGTATACCTAATGTTAAATGACGAGTTACTGTGTGCAGCAGACCAACATGGCACATGTATACATATGTAACTAACCTGCACATTGTGCGCATGTACCCGAAAACGTAAAGTATAATAAAAAAAGAAGAAGAAGAGAAAAAAAAAAGAAAAGTAAGGTGGATTAATGGATAGAGAAATGTGATAAACTGTCATAGCTAAATAATAACAATTGTAAAGCCTTGGTGATGTATATATTGGTTTTTGCTAGTATCTCTAAAAAATTTCTGGTCGAAATTTGCATAATATTCTGTGTGAAAAATCTTCTCAGGGAACATAGACCTGTGTTATTGTATTTATGTAAGTAGTTGTCATATAATTCTTTTATCTTAGAGAATCTGAAGTTTCCAATGTGAGTGTACAGTGTCTCACCACCACCTGAGTTTATTTCTGGTTGTTGTAATGCAGTACAATTGATAATGAAATTTTTTTGATGTTTTTGCCGTAGACATGGAATGTATAGTGCTCTAGCTTTTACTTTTTGAGAGACAGAAAAAGGCCCTGAGCATTGGTGTTGTGTCAGTCATCCTAGGCTCAGATCTCCATTCTACCTCTCACTGAATGCTTGGCCTTAAGCAGGTATTTAACTTTTTTTTAGATCCAATTGGTTTGTCTGTTAAATGCCTTCAGTAACGTATTCCTCACATAGTTGCTTTGAAGATTGGTGAGTCATGGCACAGCACACAGCTGAACATAGAAGACATTCAAGCAATGCCAGATTCCTTTAGGGTGCAGGTAGACATTCCTAGGGTTTGTAATAGATTTAGAAGTAGAATTAGTCTTCTACTTCTTCAGATCTCCTTTGAGGTCTGAGTGACTACAGAACCCATGGACACCTTTATTTACGAAATGCTGTGTTTCTGGTGTCTAAGGCCTGATTGCATCTTGATTGGAAACTATTAGTTAAACATTTATAAGAAATAATACTGTAACTTTTGTTTAAGCTAAAAATTATTTTCTTGAATTTTAGTTTATTGGTCTTAGTTAATACTATAACCTTTACAGTGCTAAAAAGTTATAATCCAATTGTGTTCTCTGGTGACATTAACTGAAATACATGTTCTTTTCTAGAGACCTTATTCTGTAAATCATAATAAAAAATTTCTTATGAAATTAAGTCTTTATAGATACCTCATTTAATTCTTACCATTCTATTAATAAAATACATTTTTATATTACTTATATGTATTGAGTATCTTGTCTGAAATCACACAACTAGAAAGTAGAAAGTTGTCATTTGAACACAACTCTGTGTGATTCTTAACCTTGTTTTCTTAGCCATTTTGCTGAATGGCCTTTCCTTCTTTCCTCTATCCCCCATCCCCCATCCCCCATTTCCATGTATGTGTCTGTTGTCATTCTCTACTAGCTTTTTCTTCTTATGTACAAATAGGTCAGATGTCTACATTTAAATAAACAGCATCTTTAACCCTGCTGCTCCCATACTCAAGACAAATTGACAAGTATGTTATTATCTGTGCTGTGAAGTGAAATTAGATCTTTTGATTTCAACATCAAAACCTTATAAGATGTTTTATCTTTTTTAGAAGTCCTTAATCACAAGTTAGAGCAATTAACACCTATATTACTAATAATATTTCCAGTAGTATGGATGCTTTGGTTATGGTTAAATGGTTCTGAAATCTACACCTGAGATGCTTCCCTTAAAGGGGCTCTTTGGCCAGTTCACCAATTACTTGGAAGAGGAAGAGGAAGTCTAGGTTGATTATATTGCTTTGAGTTAGAATATCTGAGAAAGTGGTACTTATATTAGTTTGCTAGTGCTGCGATGATAAAGTACTACATACTGGATGGCTTCAATAACAGAACTTTATTTTCTCACAGTCCTAGAGGCTGGAAGATCGTCTGGGACAAGTCCAGGTGCTGGCAGGGTTGATTTCTCCTGAGGCCTCTCTCTGTGGCTTGCAAGGTGGCTGTGTCTTGACATGGCCTCTTCTCTGTGTGTCCCTCCCCCTTCTTATAAGGATGTCAGTCAGATTGGATTAGGGCCCACCCTAACGGCAGCATTTTAACATTCACTTATTTAAAGACCTTATCTCTAAGAACATTACATTCTGAGGACCTGAGACTTTCAGTGTATAAATTTGGGAGGACACAATTTAGCCTATAATAGTGCTATAAAGAGAAGTTCAAAAATGAAGATGGAGATTTTATTTTTTAAGCAAAAAAATAAGTTTG";
    int M = query.length(), N = reference.length();
    auto matches = find_best_matches_np(query, reference);
    get_answer_np(matches);
    get_detail_np(matches,M,N);

    return 0;
}


