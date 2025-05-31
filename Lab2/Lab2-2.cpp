#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <functional> // For std::hash

#include <fstream>
using namespace std;
// Hash function for std::pair<int, int>
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// Define segment structure equivalent to seg_dtype
struct Segment {
    int32_t x0;       // Starting x coordinate
    int32_t y0;       // Starting y coordinate
    int32_t length;   // Segment length
    int8_t direction; // 1 for main diagonal, 0 for anti-diagonal
    int32_t distance; // Number of errors
};

// Global variables to mimic Python's global grid and tables
std::vector<std::vector<uint8_t>> grid;
std::vector<std::vector<int16_t>> main_ld, anti_lu;

// Check if two bases are complementary
bool is_complement(char a, char b) {
    return (a == 'A' && b == 'T') || (a == 'T' && b == 'A') ||
           (a == 'C' && b == 'G') || (a == 'G' && b == 'C');
}

// 构建点阵图，比较查询序列和参考序列
std::vector<std::vector<uint8_t>> build_dotplot(const std::string& query, const std::string& reference) {
    int M = query.length(), N = reference.length();
    // 初始化 M x N 的点阵图，初始值为 0
    std::vector<std::vector<uint8_t>> result(M, std::vector<uint8_t>(N, 0));
    
    for (int x = 0; x < M; ++x) {
        for (int y = 0; y < N; ++y) {
            if (query[x] == reference[y]) {
                result[x][y] = 1; // 碱基相同为1
            } else if (is_complement(query[x], reference[y])) {
                result[x][y] = 2; // 碱基互补为0
            }
        }
    }
    return result;
}

// Initialize diagonal tables
std::pair<std::vector<std::vector<int16_t>>, std::vector<std::vector<int16_t>>>
init_diag_tables(const std::vector<std::vector<uint8_t>>& grid) {
    int M = grid.size(), N = grid[0].size();
    std::vector<std::vector<int16_t>> main_ld(M, std::vector<int16_t>(N, 0));
    std::vector<std::vector<int16_t>> anti_lu(M, std::vector<int16_t>(N, 0));

    // Main diagonal: bottom-left to top-right
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            int16_t cnt = (grid[i][j] == 1) ? 1 : 0;
            if (i > 0 && j > 0) cnt += main_ld[i-1][j-1];
            main_ld[i][j] = cnt;
        }
    }

    // Anti-diagonal: top-right to bottom-left
    for (int i = 0; i < M; ++i) {
        for (int j = N-1; j >= 0; --j) {
            int16_t cnt = (grid[i][j] == 2) ? 1 : 0;
            if (i > 0 && j + 1 < N) cnt += anti_lu[i-1][j+1];
            anti_lu[i][j] = cnt;
        }
    }
    return {main_ld, anti_lu};
}

// Find diagonal segments
std::vector<Segment> find_diagonal_segments_np(const std::vector<std::vector<uint8_t>>& grid) {
    int M = grid.size(), N = grid[0].size();
    std::vector<Segment> segs;

    // Helper function to extract runs from a diagonal
    auto extract_runs = [](const std::vector<uint8_t>& line) {
        std::vector<int> starts, ends;
        std::vector<int8_t> mask(line.size() + 2, 0);
        for (size_t i = 0; i < line.size(); ++i) mask[i + 1] = (line[i] != 0);
        for (size_t i = 1; i < mask.size() - 1; ++i) {
            if (mask[i] == 1 && mask[i-1] == 0) starts.push_back(i-1);
            if (mask[i] == 1 && mask[i+1] == 0) ends.push_back(i);
        }
        return std::make_pair(starts, ends);
    };

    // Main diagonals
    for (int d = -M + 1; d < N; ++d) {
        std::vector<uint8_t> diag;
        int start_x = std::max(0, -d), start_y = start_x + d;
        for (int x = start_x, y = start_y; x < M && y < N; ++x, ++y) {
            diag.push_back(grid[x][y]);
        }
        auto [starts, ends] = extract_runs(diag);
        for (size_t i = 0; i < starts.size(); ++i) {
            int length = ends[i] - starts[i];
            if (length > 0) {
                int x0 = start_x + starts[i];
                int y0 = x0 + d;
                segs.push_back({x0, y0, length, 1, 0});
            }
        }
    }

    // Anti-diagonals (flip grid horizontally)
    std::vector<std::vector<uint8_t>> flipped(M, std::vector<uint8_t>(N));
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            flipped[i][j] = grid[i][N-1-j];
        }
    }
    for (int d = -M + 1; d < N; ++d) {
        std::vector<uint8_t> diag;
        int start_x = std::max(0, -d), start_y = N - 1 - (start_x + d);
        for (int x = start_x, y = start_y; x < M && y >= 0; ++x, --y) {
            diag.push_back(flipped[x][N-1-y]);
        }
        auto [starts, ends] = extract_runs(diag);
        for (size_t i = 0; i < starts.size(); ++i) {
            int length = ends[i] - starts[i];
            if (length > 0) {
                int x0 = start_x + starts[i];
                int y0 = N - 1 - (start_x + starts[i] + d);
                segs.push_back({x0, y0, length, 0, 0});
            }
        }
    }
    return segs;
}

// Weighted interval scheduling
std::vector<Segment> weighted_interval_scheduling_np(const std::vector<Segment>& segs) {
    if (segs.empty()) return {};

    // Sort by x1
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 + segs[a].length - 1 < segs[b].x0 + segs[b].length - 1;
    });

    std::vector<int> x0s(segs.size()), x1s(segs.size()), lens(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        x0s[i] = segs[order[i]].x0;
        lens[i] = segs[order[i]].length;
        x1s[i] = x0s[i] + lens[i] - 1;
    }

    // Compute p array
    std::vector<int> p(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        auto it = std::lower_bound(x1s.begin(), x1s.end(), x0s[i]);
        p[i] = (it == x1s.begin()) ? -1 : std::distance(x1s.begin(), it) - 1;
    }

    // Dynamic programming
    std::vector<int> dp(segs.size() + 1, 0);
    std::vector<bool> choose(segs.size(), false);
    for (size_t j = 1; j <= segs.size(); ++j) {
        int w = lens[j-1];
        if (w + dp[p[j-1] + 1] > dp[j-1]) {
            dp[j] = w + dp[p[j-1] + 1];
            choose[j-1] = true;
        } else {
            dp[j] = dp[j-1];
        }
    }

    // Backtrack
    std::vector<Segment> result;
    int j = segs.size();
    while (j > 0) {
        if (choose[j-1]) {
            result.push_back(segs[order[j-1]]);
            j = p[j-1] + 1;
        } else {
            --j;
        }
    }
    std::reverse(result.begin(), result.end());
    return result;
}

// Merge segments in blanks
std::vector<Segment> merge_in_blanks_np(const std::vector<Segment>& segs, float rate) {
    if (segs.empty()) return {};

    // Sort by x0
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 < segs[b].x0;
    });

    std::vector<Segment> merged;
    auto prev = segs[order[0]];
    for (size_t i = 1; i < segs.size(); ++i) {
        auto s = segs[order[i]];
        if (s.direction == prev.direction) {
            bool same_diag = (s.direction == 1 && (s.y0 - s.x0) == (prev.y0 - prev.x0)) ||
                             (s.direction == 0 && (s.x0 + s.y0) == (prev.x0 + prev.y0));
            if (same_diag) {
                int new_len = (s.x0 + s.length - 1) - prev.x0 + 1;
                int prev_count, affective;
                if (prev.direction == 1) {
                    prev_count = (prev.x0 > 0 && prev.y0 > 0) ? main_ld[prev.x0-1][prev.y0-1] : 0;
                    int sx1 = s.x0 + s.length - 1;
                    int sy1 = s.y0 + (s.length - 1);
                    affective = main_ld[sx1][sy1] - prev_count;
                } else {
                    prev_count = (prev.x0 > 0 && prev.y0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[prev.x0-1][prev.y0+1] : 0;
                    int sx1 = s.x0 + s.length - 1;
                    int sy1 = s.y0 - (s.length - 1);
                    affective = anti_lu[sx1][sy1] - prev_count;
                }
                int distance = new_len - affective;
                if (distance / static_cast<float>(new_len) < rate) {
                    prev.length = new_len;
                    prev.distance = distance;
                    continue;
                }
            }
        }
        merged.push_back(prev);
        prev = s;
    }
    merged.push_back(prev);
    return merged;
}

// Find large valid segments in range
std::vector<Segment> find_large_valid_segments_in_range_np(int x_start, int x_end, float rate, int min_len) {
    int M = grid.size(), N = grid[0].size();
    std::vector<Segment> recs;
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int>>, pair_hash> groups;

    // Main diagonals
    for (int d = -(M-1); d < N; ++d) {
        int lo = std::max({0, -d, x_start});
        int hi = std::min({M-1, N-1-d, x_end});
        if (hi - lo + 1 >= 1) {
            std::vector<std::pair<int, int>> pts;
            for (int x = lo; x <= hi; ++x) {
                pts.emplace_back(x, x + d);
            }
            groups[std::make_pair(1, d)] = pts;
        }
    }

    // Anti-diagonals
    for (int d = 0; d < M + N - 1; ++d) {
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

    // Process each diagonal
    for (const auto& group : groups) {
        int direction = group.first.first;
        const auto& pts = group.second;
        int L = pts.size();
        std::vector<int> mis(L + 1, 0);
        for (int i = 0; i < L; ++i) {
            auto [x, y] = pts[i];
            mis[i+1] = mis[i] + ((direction == 1 && grid[x][y] != 1) || (direction == 0 && grid[x][y] != 2));
        }

        std::vector<int> dp(L + 1, 0);
        std::vector<int> prev_idx(L + 1, 0);
        int L_ptr = 0;
        for (int R = 1; R <= L; ++R) {
            while (L_ptr < R && (mis[R] - mis[L_ptr]) > rate * (R - L_ptr)) {
                ++L_ptr;
            }
            dp[R] = dp[R-1];
            prev_idx[R] = R-1;
            if (R - L_ptr >= min_len) {
                int length = R - L_ptr;
                if (dp[L_ptr] + length > dp[R]) {
                    dp[R] = dp[L_ptr] + length;
                    prev_idx[R] = L_ptr;
                }
            }
        }

        int idx = L;
        while (idx > 0) {
            int pi = prev_idx[idx];
            if (pi < idx - 1) {
                int start = pi, end = idx - 1;
                auto [x0, y0] = pts[start];
                int length = end - start + 1;
                int errs = mis[end + 1] - mis[start];
                recs.push_back(Segment{x0, y0, length, static_cast<int8_t>(direction), errs});
                idx = pi;
            } else {
                --idx;
            }
        }
    }
    return recs;
}

// Fill in blanks globally
std::vector<Segment> fill_in_blanks_global_np(const std::vector<Segment>& segs, float rate, int min_gap) {
    if (segs.empty()) return {};

    // Sort by x0
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 < segs[b].x0;
    });

    std::vector<Segment> recs;
    for (size_t i = 0; i < segs.size() - 1; ++i) {
        auto prev = segs[order[i]];
        auto curr = segs[order[i + 1]];
        recs.push_back(prev);
        int g0 = prev.x0 + prev.length;
        int g1 = curr.x0 - 1;
        int gap_len = g1 - g0 + 1;
        if (gap_len >= min_gap) {
            auto extras = find_large_valid_segments_in_range_np(g0, g1, rate, min_gap);
            recs.insert(recs.end(), extras.begin(), extras.end());
        }
    }
    recs.push_back(segs[order.back()]);

    recs = merge_in_blanks_np(recs, rate);
    recs = weighted_interval_scheduling_np(recs);
    return recs;
}

// Extend segment ends backward
std::vector<Segment> extend_end_backward_np(const std::vector<Segment>& segs, float rate) {
    if (segs.empty()) return {};

    // Sort by x0 descending
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 > segs[b].x0;
    });

    std::vector<Segment> out;
    int prev_end = segs[order[0]].x0 + segs[order[0]].length;
    for (size_t i = 0; i < segs.size(); ++i) {
        auto seg = segs[order[i]];
        int ox0 = seg.x0, oy0 = seg.y0, olen = seg.length, odir = seg.direction, odist = seg.distance;
        int ox1 = ox0 + olen - 1;
        int oy1 = oy0 + (olen - 1) * (odir == 1 ? 1 : -1);

        int target_end = prev_end - 1;
        int space = target_end - ox1;
        int new_len = olen, new_dist = odist;
        int cand_x1 = ox1, cand_y1 = oy1;

        if (ox1 < target_end) {
            cand_x1 = target_end;
            int step = (odir == 1) ? 1 : -1;
            cand_y1 = oy1 + step * space;
            new_len = olen + space;

            while (cand_x1 > ox1) {
                if (cand_x1 >= static_cast<int>(grid.size()) || cand_y1 < 0 || cand_y1 >= static_cast<int>(grid[0].size())) {
                    --cand_x1;
                    cand_y1 -= step;
                    --new_len;
                    continue;
                }
                int prev_corr, corr;
                if (odir == 1) {
                    prev_corr = (ox0 > 0 && oy0 > 0) ? main_ld[ox0-1][oy0-1] : 0;
                    corr = main_ld[cand_x1][cand_y1] - prev_corr;
                } else {
                    prev_corr = (ox0 > 0 && oy0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[ox0-1][oy0+1] : 0;
                    corr = anti_lu[cand_x1][cand_y1] - prev_corr;
                }
                int dist = new_len - corr;
                if (dist / static_cast<float>(new_len) < rate) {
                    new_dist = dist;
                    break;
                }
                --cand_x1;
                cand_y1 -= step;
                --new_len;
            }
        }
        out.push_back(Segment{ox0, oy0, new_len, static_cast<int8_t>(odir), new_dist});
        prev_end = ox0;
    }

    // Sort back by x0 ascending
    std::vector<size_t> asc_order(out.size());
    for (size_t i = 0; i < out.size(); ++i) asc_order[i] = i;
    std::sort(asc_order.begin(), asc_order.end(), [&out](size_t a, size_t b) {
        return out[a].x0 < out[b].x0;
    });
    std::vector<Segment> result(out.size());
    for (size_t i = 0; i < out.size(); ++i) result[i] = out[asc_order[i]];
    return result;
}

// Extend segment starts backward
std::vector<Segment> extend_start_backward_np(const std::vector<Segment>& segs, float rate) {
    if (segs.empty()) return {};

    // Sort by x0 ascending
    std::vector<size_t> order(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&segs](size_t a, size_t b) {
        return segs[a].x0 < segs[b].x0;
    });

    std::vector<Segment> out;
    int prev_end = -1;
    for (size_t i = 0; i < segs.size(); ++i) {
        auto seg = segs[order[i]];
        int ox0 = seg.x0, oy0 = seg.y0, olen = seg.length, odir = seg.direction, odist = seg.distance;
        int ox1 = ox0 + olen - 1;
        int oy1 = oy0 + (olen - 1) * (odir == 1 ? 1 : -1);

        int target_start = prev_end + 1;
        int space = ox0 - target_start;
        int cand_x0 = ox0, cand_y0 = oy0, cand_len = olen;

        if (ox0 > target_start) {
            cand_x0 = target_start;
            int step = (odir == 1) ? 1 : -1;
            cand_y0 = oy0 - step * space;
            cand_len = olen + space;

            while (cand_x0 < ox0) {
                if (cand_x0 < 0 || cand_y0 < 0 || cand_y0 >= static_cast<int>(grid[0].size())) break;
                int prev_corr, corr;
                if (odir == 1) {
                    prev_corr = (cand_x0 > 0 && cand_y0 > 0) ? main_ld[cand_x0-1][cand_y0-1] : 0;
                    corr = main_ld[cand_x0 + cand_len - 1][cand_y0 + cand_len - 1] - prev_corr;
                } else {
                    prev_corr = (cand_x0 > 0 && cand_y0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[cand_x0-1][cand_y0+1] : 0;
                    corr = anti_lu[cand_x0 + cand_len - 1][cand_y0 - (cand_len - 1)] - prev_corr;
                }
                int dist = cand_len - corr;
                if (dist / static_cast<float>(cand_len) < rate) {
                    out.push_back(Segment{cand_x0, cand_y0, cand_len, static_cast<int8_t>(odir), dist});
                    break;
                }
                ++cand_x0;
                cand_y0 += step;
                --cand_len;
            }
            if (cand_x0 == ox0) out.push_back(seg);
        } else {
            out.push_back(seg);
        }
        prev_end = out.back().x0 + out.back().length - 1;
    }
    return out;
}

// Choose segments with length >= min_len
std::vector<Segment> chose_segs_np(const std::vector<Segment>& segs, int min_len) {
    std::vector<Segment> result;
    for (const auto& seg : segs) {
        if (seg.length >= min_len) result.push_back(seg);
    }
    return result;
}

// Minimal interval cover
std::vector<Segment> minimal_interval_cover2_np(const std::vector<Segment>& segs, float rate, int length_thresh) {
    if (segs.empty()) return {};

    // Find global min_x and max_x
    int min_x = segs[0].x0, max_x = segs[0].x0 + segs[0].length - 1;
    for (const auto& s : segs) {
        min_x = std::min(min_x, s.x0);
        max_x = std::max(max_x, s.x0 + s.length - 1);
    }

    // Build buckets
    std::vector<int> best_end(max_x + 1, -1);
    std::vector<int> best_idx(max_x + 1, -1);
    for (size_t idx = 0; idx < segs.size(); ++idx) {
        int bx0 = segs[idx].x0;
        int bx1 = bx0 + segs[idx].length - 1;
        if (bx1 > best_end[bx0]) {
            best_end[bx0] = bx1;
            best_idx[bx0] = idx;
        }
    }

    std::vector<Segment> result;
    int covered_end = min_x - 1;
    while (covered_end < max_x) {
        int start_pos = covered_end + 1;
        int candidate_end = -1, candidate_i = -1;
        for (int x0 = min_x; x0 <= start_pos; ++x0) {
            if (best_end[x0] > candidate_end) {
                candidate_end = best_end[x0];
                candidate_i = best_idx[x0];
            }
        }

        if (candidate_i < 0 || candidate_end < start_pos) break;

        auto best = segs[candidate_i];
        int new_len = std::min(best.x0 + best.length - 1, candidate_end) - start_pos + 1;
        int new_y0, prev, correct;
        if (best.direction == 1) {
            new_y0 = best.y0 + (start_pos - best.x0);
            prev = (start_pos > 0 && new_y0 > 0) ? main_ld[start_pos-1][new_y0-1] : 0;
            correct = main_ld[start_pos + new_len - 1][new_y0 + new_len - 1] - prev;
        } else {
            new_y0 = best.y0 - (start_pos - best.x0);
            prev = (start_pos > 0 && new_y0 + 1 < static_cast<int>(grid[0].size())) ? anti_lu[start_pos-1][new_y0+1] : 0;
            correct = anti_lu[start_pos + new_len - 1][new_y0 - new_len + 1] - prev;
        }
        int distance = new_len - correct;

        if (new_len >= length_thresh && distance / static_cast<float>(new_len) < rate) {
            result.push_back(Segment{start_pos, new_y0, new_len, static_cast<int8_t>(best.direction), distance});
        }
        covered_end = candidate_end;
    }
    return result;
}

// Merge segments with tolerance
std::vector<Segment> merge_with_tolerance_np(const std::vector<Segment>& segs, int max_gap, float max_error_rate) {
    if (segs.empty()) return {};

    // Group by direction and diagonal id
    std::unordered_map<std::pair<int, int>, std::vector<size_t>, pair_hash> groups;
    std::vector<int> diag_id(segs.size());
    std::vector<int> x1(segs.size()), y1(segs.size());

    for (size_t i = 0; i < segs.size(); ++i) {
        x1[i] = segs[i].x0 + segs[i].length - 1;
        y1[i] = segs[i].direction == 1 ? segs[i].y0 + segs[i].length - 1 : segs[i].y0 - (segs[i].length - 1);
        diag_id[i] = segs[i].direction == 1 ? segs[i].y0 - segs[i].x0 : segs[i].y0 + segs[i].x0;
        groups[std::make_pair(segs[i].direction, diag_id[i])].push_back(i);
    }

    std::vector<Segment> merged;
    for (const auto& group : groups) {
        int dir_val = group.first.first;
        auto indices = group.second;

        // Sort group by x0
        std::sort(indices.begin(), indices.end(), [&segs](size_t a, size_t b) {
            return segs[a].x0 < segs[b].x0;
        });

        auto curr = segs[indices[0]];
        int cx1 = curr.x0 + curr.length - 1;
        int cy1 = curr.direction == 1 ? curr.y0 + curr.length - 1 : curr.y0 - (curr.length - 1);

        for (size_t i = 1; i < indices.size(); ++i) {
            auto s = segs[indices[i]];
            int sx1 = s.x0 + s.length - 1;
            int sy1 = s.direction == 1 ? s.y0 + s.length - 1 : s.y0 - (s.length - 1);
            int gap = s.x0 - cx1 - 1;
            int merged_len = sx1 - curr.x0 + 1;
            bool cond_align = (dir_val == 1 && s.y0 - cy1 == gap + 1) ||
                              (dir_val == 0 && cy1 - s.y0 == gap + 1);

            if (gap <= max_gap && (gap + curr.distance + s.distance) / static_cast<float>(merged_len) <= max_error_rate && cond_align) {
                curr.length = merged_len;
                curr.distance = gap + curr.distance + s.distance;
                cx1 = curr.x0 + curr.length - 1;
                cy1 = curr.y0 + (curr.length - 1) * (curr.direction == 1 ? 1 : -1);
            } else {
                merged.push_back(curr);
                curr = s;
                cx1 = curr.x0 + curr.length - 1;
                cy1 = curr.y0 + (curr.length - 1) * (curr.direction == 1 ? 1 : -1);
            }
        }
        merged.push_back(curr);
    }
    return merged;
}

// Print answer
void get_answer_np(const std::vector<Segment>& matches) {
    std::cout << "Remain segments:" << std::endl;
    for (const auto& seg : matches) {
        if (seg.length >= 30) {
            int x0 = seg.x0, length = seg.length, y0 = seg.y0;
            int x1 = x0 + length - 1;
            if (seg.direction == 1) {
                int y1 = y0 + length - 1;
                std::cout << "(" << x0 << "," << x1 + 1 << "," << y0 << "," << y1 + 1 << "),";
            } else {
                int y1 = y0 - (length - 1);
                std::cout << "(" << x0 << "," << x1 + 1 << "," << y1 << "," << y0 + 1 << "),";
            }
        }
    }
    std::cout << std::endl;
}

// Print detailed output
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

// Main function
std::vector<Segment> find_best_matches_np(const std::string& query, const std::string& reference) {
    grid = build_dotplot(query, reference); // 构建点阵图，比较查询序列和参考序列
    auto [ml, al] = init_diag_tables(grid); // 初始化主对角线和反对角线匹配表
    main_ld = ml;
    anti_lu = al;

    auto final_segs = find_diagonal_segments_np(grid); // 提取对角线上的匹配段
    final_segs = merge_with_tolerance_np(final_segs, 1, 0.039); // 合并具有小间隙的段
    final_segs = minimal_interval_cover2_np(final_segs, 0.08, 20); // 最小区间覆盖，选择最少段覆盖整个区间
    final_segs = fill_in_blanks_global_np(final_segs, 0.065, 25); // 全局填充空白区域
    final_segs = chose_segs_np(final_segs, 25); // 选择长度大于等于 min_len 的段
    final_segs = extend_start_backward_np(final_segs, 0.1); // 向前扩展段的起始位置

    // Sort by x0 descending for extend_end_backward
    std::vector<size_t> order(final_segs.size());
    for (size_t i = 0; i < final_segs.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&final_segs](size_t a, size_t b) {
        return final_segs[a].x0 > final_segs[b].x0;
    });
    std::vector<Segment> temp(final_segs.size());
    for (size_t i = 0; i < final_segs.size(); ++i) temp[i] = final_segs[order[i]];
    final_segs = extend_end_backward_np(temp, 0.1);  // 向后扩展段的结束位置

    // Sort back by x0 ascending
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
