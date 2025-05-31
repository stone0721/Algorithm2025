#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <functional>

// Hash function for std::pair<int, int>
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// Define segment structure with int16_t for coordinates
struct Segment {
    int16_t x0;       // Starting x coordinate
    int16_t y0;       // Starting y coordinate
    int16_t length;   // Segment length
    int8_t direction; // 1 for main diagonal, 0 for anti-diagonal
    int16_t distance; // Number of errors
};

// Global sequences (avoid copying)
const std::string* query_ptr = nullptr;
const std::string* reference_ptr = nullptr;

// Compute grid value on-the-fly
uint8_t get_grid_value(int x, int y) {
    char q = (*query_ptr)[x], r = (*reference_ptr)[y];
    if (q == r) return 1;
    if ((q == 'A' && r == 'T') || (q == 'T' && r == 'A') ||
        (q == 'C' && r == 'G') || (q == 'G' && r == 'C')) return 2;
    return 0;
}

// Compute matches in a diagonal range
int compute_diagonal_matches(int x_start, int y_start, int length, int direction, int M, int N) {
    int matches = 0;
    for (int i = 0; i < length; ++i) {
        int x = x_start + i;
        int y = direction == 1 ? y_start + i : y_start - i;
        if (x >= M || y < 0 || y >= N) break;
        uint8_t val = get_grid_value(x, y);
        matches += (direction == 1 && val == 1) || (direction == 0 && val == 2);
    }
    return matches;
}

// Find diagonal segments without storing grid
std::vector<Segment> find_diagonal_segments_np(int M, int N, int min_len = 5) {
    std::vector<Segment> segs;

    // Main diagonals
    for (int d = -M + 1; d < N; ++d) {
        int start_x = std::max(0, -d), start_y = start_x + d;
        std::vector<uint8_t> diag;
        for (int x = start_x, y = start_y; x < M && y < N; ++x, ++y) {
            diag.push_back(get_grid_value(x, y));
        }
        // Extract runs
        std::vector<int> starts, ends;
        std::vector<int8_t> mask(diag.size() + 2, 0);
        for (size_t i = 0; i < diag.size(); ++i) mask[i + 1] = (diag[i] == 1);
        for (size_t i = 1; i < mask.size() - 1; ++i) {
            if (mask[i] == 1 && mask[i-1] == 0) starts.push_back(i-1);
            if (mask[i] == 1 && mask[i+1] == 0) ends.push_back(i);
        }
        for (size_t i = 0; i < starts.size(); ++i) {
            int length = ends[i] - starts[i];
            if (length >= min_len) {
                int x0 = start_x + starts[i];
                int y0 = x0 + d;
                segs.push_back({static_cast<int16_t>(x0), static_cast<int16_t>(y0), 
                               static_cast<int16_t>(length), 1, 0});
            }
        }
    }

    // Anti-diagonals
    for (int d = -M + 1; d < N; ++d) {
        int start_x = std::max(0, -d), start_y = N - 1 - (start_x + d);
        std::vector<uint8_t> diag;
        for (int x = start_x, y = start_y; x < M && y >= 0; ++x, --y) {
            diag.push_back(get_grid_value(x, N-1-y));
        }
        // Extract runs
        std::vector<int> starts, ends;
        std::vector<int8_t> mask(diag.size() + 2, 0);
        for (size_t i = 0; i < diag.size(); ++i) mask[i + 1] = (diag[i] == 2);
        for (size_t i = 1; i < mask.size() - 1; ++i) {
            if (mask[i] == 1 && mask[i-1] == 0) starts.push_back(i-1);
            if (mask[i] == 1 && mask[i+1] == 0) ends.push_back(i);
        }
        for (size_t i = 0; i < starts.size(); ++i) {
            int length = ends[i] - starts[i];
            if (length >= min_len) {
                int x0 = start_x + starts[i];
                int y0 = N - 1 - (start_x + starts[i] + d);
                segs.push_back({static_cast<int16_t>(x0), static_cast<int16_t>(y0), 
                               static_cast<int16_t>(length), 0, 0});
            }
        }
    }
    return segs;
}

// Weighted interval scheduling
std::vector<Segment> weighted_interval_scheduling_np(const std::vector<Segment>& segs) {
    if (segs.empty()) return {};

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

    std::vector<int> p(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        auto it = std::lower_bound(x1s.begin(), x1s.end(), x0s[i]);
        p[i] = (it == x1s.begin()) ? -1 : std::distance(x1s.begin(), it) - 1;
    }

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
std::vector<Segment> merge_in_blanks_np(const std::vector<Segment>& segs, float rate, int M, int N) {
    if (segs.empty()) return {};

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
                int affective = compute_diagonal_matches(prev.x0, prev.y0, new_len, prev.direction, M, N);
                int distance = new_len - affective;
                if (distance / static_cast<float>(new_len) < rate) {
                    prev.length = static_cast<int16_t>(new_len);
                    prev.distance = static_cast<int16_t>(distance);
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
std::vector<Segment> find_large_valid_segments_in_range_np(int x_start, int x_end, float rate, int min_len, int M, int N) {
    std::vector<Segment> recs;
    std::unordered_map<std::pair<int, int>, std::vector<std::pair<int, int>>, pair_hash> groups;

    // Main diagonals
    for (int d = -(M-1); d < N; ++d) {
        int lo = std::max({0, -d, x_start});
        int hi = std::min({M-1, N-1-d, x_end});
        if (hi - lo + 1 >= min_len) {
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
        if (hi - lo + 1 >= min_len) {
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
            mis[i+1] = mis[i] + ((direction == 1 && get_grid_value(x, y) != 1) || (direction == 0 && get_grid_value(x, y) != 2));
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

// Fill in blanks globally
std::vector<Segment> fill_in_blanks_global_np(const std::vector<Segment>& segs, float rate, int min_gap, int M, int N) {
    if (segs.empty()) return {};

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
            auto extras = find_large_valid_segments_in_range_np(g0, g1, rate, min_gap, M, N);
            recs.insert(recs.end(), extras.begin(), extras.end());
        }
    }
    recs.push_back(segs[order.back()]);

    recs = merge_in_blanks_np(recs, rate, M, N);
    recs = weighted_interval_scheduling_np(recs);
    return recs;
}

// Extend segment ends backward
std::vector<Segment> extend_end_backward_np(const std::vector<Segment>& segs, float rate, int M, int N) {
    if (segs.empty()) return {};

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
                if (cand_x1 >= M || cand_y1 < 0 || cand_y1 >= N) {
                    --cand_x1;
                    cand_y1 -= step;
                    --new_len;
                    continue;
                }
                int corr = compute_diagonal_matches(ox0, oy0, new_len, odir, M, N);
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
        out.push_back({static_cast<int16_t>(ox0), static_cast<int16_t>(oy0), 
                      static_cast<int16_t>(new_len), static_cast<int8_t>(odir), 
                      static_cast<int16_t>(new_dist)});
        prev_end = ox0;
    }

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
std::vector<Segment> extend_start_backward_np(const std::vector<Segment>& segs, float rate, int M, int N) {
    if (segs.empty()) return {};

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
                if (cand_x0 < 0 || cand_y0 < 0 || cand_y0 >= N) break;
                int corr = compute_diagonal_matches(cand_x0, cand_y0, cand_len, odir, M, N);
                int dist = cand_len - corr;
                if (dist / static_cast<float>(cand_len) < rate) {
                    out.push_back({static_cast<int16_t>(cand_x0), static_cast<int16_t>(cand_y0), 
                                  static_cast<int16_t>(cand_len), static_cast<int8_t>(odir), 
                                  static_cast<int16_t>(dist)});
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
std::vector<Segment> minimal_interval_cover2_np(const std::vector<Segment>& segs, float rate, int length_thresh, int M, int N) {
    if (segs.empty()) return {};

    int min_x = segs[0].x0, max_x = segs[0].x0 + segs[0].length - 1;
    for (const auto& s : segs) {
        min_x = std::min(min_x, static_cast<int>(s.x0));
        max_x = std::max(max_x, static_cast<int>(s.x0 + s.length - 1));
    }

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
        int new_len = std::min(static_cast<int>(best.x0 + best.length - 1), candidate_end) - start_pos + 1;
        int new_y0 = best.direction == 1 ? best.y0 + (start_pos - best.x0) : best.y0 - (start_pos - best.x0);
        int correct = compute_diagonal_matches(start_pos, new_y0, new_len, best.direction, M, N);
        int distance = new_len - correct;

        if (new_len >= length_thresh && distance / static_cast<float>(new_len) < rate) {
            result.push_back({static_cast<int16_t>(start_pos), static_cast<int16_t>(new_y0), 
                            static_cast<int16_t>(new_len), static_cast<int8_t>(best.direction), 
                            static_cast<int16_t>(distance)});
        }
        covered_end = candidate_end;
    }
    return result;
}

// Merge segments with tolerance
std::vector<Segment> merge_with_tolerance_np(const std::vector<Segment>& segs, int max_gap, float max_error_rate, int M, int N) {
    if (segs.empty()) return {};

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
                curr.length = static_cast<int16_t>(merged_len);
                curr.distance = static_cast<int16_t>(gap + curr.distance + s.distance);
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
void get_detail_np(const std::vector<Segment>& matches, int M, int N) {
    int point = 0;
    for (const auto& seg : matches) {
        if (seg.length <= 0) continue;
        int x0 = seg.x0, y0 = seg.y0, length = seg.length, direction = seg.direction;
        int x1 = x0 + length - 1;
        int affective = compute_diagonal_matches(x0, y0, length, direction, M, N);
        //if (direction == 1) {
        //    int y1 = y0 + length - 1;
        //    std::cout << "(" << x0 << "," << x1 + 1 << "," << y0 << "," << y1 + 1 << ")\t";
        //} else {
        //    int y1 = y0 - (length - 1);
        //    std::cout << "(" << x0 << "," << x1 + 1 << "," << y1 << "," << y0 + 1 << ")\t";
        //}
        //std::cout << ";length = " << length << "\taffective length = " << affective << "\t";
        if (affective / static_cast<float>(length) > 0.9) {
            //std::cout << "True" << std::endl;
            point += affective;
        } else {
            //std::cout << "False" << std::endl;
        }
    }
    std::cout << "final score = " << point << std::endl;
}

// Main function
std::vector<Segment> find_best_matches_np(const std::string& query, const std::string& reference) {
    query_ptr = &query;
    reference_ptr = &reference;
    int M = query.length(), N = reference.length();

    auto final_segs = find_diagonal_segments_np(M, N, 5);
    final_segs = merge_with_tolerance_np(final_segs, 3, 0.08, M, N);
    final_segs = minimal_interval_cover2_np(final_segs, 0.085, 30, M, N);
    final_segs = fill_in_blanks_global_np(final_segs, 0.075, 15, M, N);
    final_segs = chose_segs_np(final_segs, 30);
    final_segs = extend_start_backward_np(final_segs, 0.085, M, N);
    final_segs = extend_end_backward_np(final_segs, 0.085, M, N);

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


