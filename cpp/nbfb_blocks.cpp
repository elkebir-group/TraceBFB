// Index-tuple search for the nested BFB distance over block sources; see
// nbfb_blocks.hpp and PREFIX_FAST_DESIGN.md.
//
// Text. X = rc(t) # g_1 # ... # g_r # rc(g_r) # ... # rc(g_1) # t, 1-based,
// with 0 as the separator. X is its own reverse complement, so the reverse
// complement of X[a,b] is X[N+1-b, N+1-a]. Positions 1..R (everything
// before the last separator) form D: a tail u makes t[1,phi] u a node
// exactly when u occurs in D. Positions n+2..R are the sources.
//
// Nodes. E(phi) = t[1,phi]; W(phi,a,b) = t[1,phi] X[a,b] with a..b inside
// D and X[a] != t[phi+1]. Several tuples may spell the same string.
//
// Query. occ(x, l, g) is the maximum over occurrences p of X[x,x+l-1] in D
// of LCE(p+l, g), with an argmax p, memoized. Every rule of a W node fixes
// the left part of the produced tail and lets the right part grow, so one
// query per cut point gives all produced tails as prefixes of X[p, ...]:
// an interval of tuples W(phi', p, b') over b'.

#include "nbfb_blocks.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace nbfb {

namespace {

using i32 = std::int32_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 NO_PARENT = 0xffffffffu;

std::vector<int> reverse_complement(const std::vector<int>& x) {
    std::vector<int> y(x.rbegin(), x.rend());
    for (int& v : y) v = -v;
    return y;
}

bool contains(const std::vector<int>& text, const std::vector<int>& pattern) {
    return std::search(text.begin(), text.end(), pattern.begin(),
                       pattern.end()) != text.end();
}

struct Record {
    i32 phi;
    i32 a;  // W only
    i32 b;  // W only
    u32 parent;
    u32 dist;
    std::uint16_t mask;
    char kind;     // 'E' or 'W'
    char op_kind;  // 'p', 'i', 'f', 't'
    i32 op_x;      // beta, or the first position of the appended block
    i32 op_y;      // gamma / c, or the last position of the appended block
};

// Next-unvisited structure over 0..size-1 (interval union-find).
struct Family {
    std::vector<i32> parent;
    explicit Family(i32 size) : parent(static_cast<std::size_t>(size) + 1) {
        for (i32 v = 0; v <= size; ++v) parent[static_cast<std::size_t>(v)] = v;
    }
    i32 next(i32 v) {
        i32 root = v;
        while (parent[static_cast<std::size_t>(root)] != root) {
            root = parent[static_cast<std::size_t>(root)];
        }
        while (parent[static_cast<std::size_t>(v)] != root) {
            const i32 up = parent[static_cast<std::size_t>(v)];
            parent[static_cast<std::size_t>(v)] = root;
            v = up;
        }
        return root;
    }
    void remove(i32 v) { parent[static_cast<std::size_t>(v)] = v + 1; }
};

struct Layer {
    Family e_family;  // E(phi), phi = 0..n
    std::unordered_map<u64, u32> w_index;  // (phi, start) -> family
    std::vector<Family> w_families;        // offsets b - start = 0..n-1
    explicit Layer(i32 n) : e_family(n + 1) {}
};

class BlockSearch {
  public:
    BlockSearch(const PreparedSources& sources, const std::vector<int>& target,
                MemoryBudget& budget, const JoinTable* joins)
        : target_(target), n_(static_cast<i32>(target.size())),
          budget_(budget), joins_(joins) {
        build_text(sources);
        build_lce();
        for (std::size_t k = 0; k < sources.forward.size(); ++k) {
            listed_.emplace(sources.forward[k], static_cast<int>(k));
        }
    }

    ~BlockSearch() { budget_.release(owned_); }

    Result solve(bool want_history) {
        Result result;
        long long goal = -1;
        if (joins_ == nullptr) {
            for (int bound = 1; bound <= n_ && goal < 0; ++bound) {
                reset(bound, 0, true);
                goal = run();
            }
        } else {
            const int cap = std::min<int>(n_, 65535 / stride_ - 1);
            for (int fb = 0; fb <= cap && goal < 0; ++fb) {
                reset(n_, fb, false);
                if (run() < 0) continue;
                for (int bound = 1; bound <= n_ && goal < 0; ++bound) {
                    reset(bound, fb, true);
                    goal = run();
                }
            }
            if (goal < 0 && cap < n_) {
                throw std::invalid_argument(
                    "target too long for the fold-back-first objective");
            }
        }
        if (goal < 0) return result;
        const Record& end = records_[static_cast<std::size_t>(goal)];
        result.reachable = true;
        result.t_cycles = t_of(end.mask);
        result.foldback_t = f_of(end.mask);
        result.total = static_cast<int>(end.dist);
        const std::vector<Operation> history = replay(static_cast<u32>(goal));
        if (want_history) result.history = history;
        return result;
    }

  private:
    // ------------------------------------------------------------ text

    void build_text(const PreparedSources& sources) {
        // The maximal sources: longest first, dropping every source that
        // occurs in a kept one or in the reverse complement of a kept one.
        std::vector<std::vector<int>> all = sources.forward;
        for (const std::vector<int>& x : all) {
            if (x.empty()) throw std::invalid_argument("source is empty");
        }
        std::stable_sort(all.begin(), all.end(),
                         [](const std::vector<int>& l, const std::vector<int>& r) {
                             return l.size() > r.size();
                         });
        std::vector<std::vector<int>> kept;
        for (const std::vector<int>& x : all) {
            bool inside = false;
            for (const std::vector<int>& g : kept) {
                if (contains(g, x) || contains(reverse_complement(g), x)) {
                    inside = true;
                    break;
                }
            }
            if (!inside) kept.push_back(x);
        }

        X_.assign(1, 0);  // X[0] unused
        const std::vector<int> tbar = reverse_complement(target_);
        X_.insert(X_.end(), tbar.begin(), tbar.end());
        for (const std::vector<int>& g : kept) {
            X_.push_back(0);
            X_.insert(X_.end(), g.begin(), g.end());
        }
        for (auto it = kept.rbegin(); it != kept.rend(); ++it) {
            X_.push_back(0);
            const std::vector<int> gbar = reverse_complement(*it);
            X_.insert(X_.end(), gbar.begin(), gbar.end());
        }
        X_.push_back(0);
        R_ = static_cast<i32>(X_.size()) - 2;  // last position of D
        X_.insert(X_.end(), target_.begin(), target_.end());
        N_ = static_cast<i32>(X_.size()) - 1;
        if (N_ >= (1 << 20)) throw std::invalid_argument("instance too large");
        stride_ = n_ + 1;
    }

    void build_lce() {
        width_ = static_cast<std::size_t>(N_) + 2;
        const std::size_t bytes = width_ * width_ * sizeof(i32);
        acquire(bytes);
        lce_.assign(width_ * width_, 0);
        for (i32 i = N_; i >= 1; --i) {
            for (i32 j = N_; j >= 1; --j) {
                const int x = X_[static_cast<std::size_t>(i)];
                if (x != 0 && x == X_[static_cast<std::size_t>(j)]) {
                    lce_[idx(i, j)] = 1 + lce_[idx(i + 1, j + 1)];
                }
            }
        }
    }

    std::size_t idx(i32 i, i32 j) const {
        return static_cast<std::size_t>(i) * width_ + static_cast<std::size_t>(j);
    }
    i32 lce(i32 i, i32 j) const { return lce_[idx(i, j)]; }
    i32 tpos(i32 k) const { return N_ - n_ + k; }  // X position of t[k]

    // Maximum over occurrences p of X[x, x+l-1] inside D of LCE(p+l, g).
    // Returns -1 if the window does not occur in D.
    i32 occ(i32 x, i32 l, i32 g, i32& pos) {
        const u64 key = (static_cast<u64>(x) << 42) | (static_cast<u64>(l) << 21) |
                        static_cast<u64>(g);
        auto found = occ_.find(key);
        if (found != occ_.end()) {
            pos = found->second.second;
            return found->second.first;
        }
        i32 best = -1;
        pos = -1;
        for (i32 p = 1; p + l - 1 <= R_; ++p) {
            if (lce(p, x) >= l) {
                const i32 value = lce(p + l, g);
                if (value > best) {
                    best = value;
                    pos = p;
                }
            }
        }
        occ_.emplace(key, std::make_pair(best, pos));
        if (occ_.size() % 4096 == 0) acquire(4096 * 48);
        return best;
    }

    // ------------------------------------------------------------ masks

    int t_of(int state) const { return state % stride_; }
    int f_of(int state) const { return state / stride_; }
    bool may_append(int state) const { return t_of(state) < max_t_; }
    int after_t(int state, int u, int v) const {
        const int fb =
            (joins_ != nullptr && u != 0 && joins_->is_foldback(u, v)) ? 1 : 0;
        if (f_of(state) + fb > max_fb_) return -1;
        return state + (count_t_ ? 1 : 0) + fb * stride_;
    }

    void reset(int max_t, int max_fb, bool count_t) {
        max_t_ = max_t;
        max_fb_ = max_fb;
        count_t_ = count_t;
        records_.clear();
        layers_.clear();
        budget_.release(search_bytes_);
        owned_ -= search_bytes_;
        search_bytes_ = 0;
    }

    Layer& layer(int mask) {
        auto found = layers_.find(mask);
        if (found != layers_.end()) return found->second;
        charge(static_cast<std::size_t>(n_ + 2) * sizeof(i32) + 64);
        return layers_.emplace(mask, Layer(n_)).first->second;
    }

    Family& w_family(Layer& L, i32 phi, i32 start) {
        const u64 key = (static_cast<u64>(phi) << 32) | static_cast<u64>(start);
        auto found = L.w_index.find(key);
        if (found != L.w_index.end()) return L.w_families[found->second];
        charge(static_cast<std::size_t>(n_ + 1) * sizeof(i32) + 64);
        L.w_index.emplace(key, static_cast<u32>(L.w_families.size()));
        L.w_families.emplace_back(n_);
        return L.w_families.back();
    }

    void push(char kind, int mask, i32 phi, i32 a, i32 b, u32 parent,
              char op_kind, i32 op_x, i32 op_y) {
        if (records_.size() == records_.capacity()) {
            charge((records_.capacity() + 1024) * sizeof(Record));
        }
        Record r;
        r.phi = phi;
        r.a = a;
        r.b = b;
        r.parent = parent;
        r.dist = parent == NO_PARENT ? 0u : records_[parent].dist + 1u;
        r.mask = static_cast<std::uint16_t>(mask);
        r.kind = kind;
        r.op_kind = op_kind;
        r.op_x = op_x;
        r.op_y = op_y;
        records_.push_back(r);
    }

    void visit_E(int mask, i32 phi, u32 parent, char op, i32 x, i32 y) {
        if (phi > n_) return;
        Family& f = layer(mask).e_family;
        if (f.next(phi) != phi) return;
        f.remove(phi);
        push('E', mask, phi, 0, 0, parent, op, x, y);
    }

    // Visits E(v) for v in [lo, hi]; op_of(v) gives (beta, cut).
    template <class OpOf>
    void visit_E_interval(int mask, i32 lo, i32 hi, u32 parent, char op,
                          OpOf op_of) {
        hi = std::min(hi, n_);
        if (lo > hi) return;
        Family& f = layer(mask).e_family;
        for (i32 v = f.next(lo); v <= hi; v = f.next(v + 1)) {
            f.remove(v);
            i32 x, y;
            op_of(v, x, y);
            push('E', mask, v, 0, 0, parent, op, x, y);
        }
    }

    // Visits W(phi, start, b) for b in [lo, hi], clipped to length <= n;
    // op_of(b) gives the cycle parameters.
    template <class OpOf>
    void visit_W_interval(int mask, i32 phi, i32 start, i32 lo, i32 hi,
                          u32 parent, char op, OpOf op_of) {
        lo = std::max(lo, start);
        hi = std::min(hi, n_ - phi + start - 1);
        if (lo > hi) return;
        Layer& L = layer(mask);
        Family& f = w_family(L, phi, start);
        for (i32 off = f.next(lo - start); off <= hi - start; off = f.next(off + 1)) {
            f.remove(off);
            i32 x, y;
            op_of(start + off, x, y);
            push('W', mask, phi, start, start + off, parent, op, x, y);
        }
    }

    // ------------------------------------------------------------ search

    long long run() {
        push('E', 0, 0, 0, 0, NO_PARENT, 'p', -1, -1);
        layer(0).e_family.remove(0);
        for (std::size_t head = 0; head < records_.size(); ++head) {
            const Record r = records_[head];  // push may reallocate
            if (r.kind == 'E') {
                if (r.phi == n_) return static_cast<long long>(head);
                expand_E(r, static_cast<u32>(head));
            } else {
                expand_W(r, static_cast<u32>(head));
            }
        }
        return -1;
    }

    void expand_E(const Record& r, u32 key) {
        const int mask = r.mask;
        const i32 phi = r.phi;
        const i32 n = n_;
        // t-BFB: append a block X[i,j] of the sources.
        if (may_append(mask)) {
            const int last = phi >= 1 ? target_[static_cast<std::size_t>(phi - 1)] : 0;
            for (i32 i = n + 2; i <= R_; ++i) {
                if (X_[static_cast<std::size_t>(i)] == 0) continue;
                const int next = after_t(mask, last, X_[static_cast<std::size_t>(i)]);
                if (next < 0) continue;
                const i32 m = lce(i, tpos(phi + 1));
                for (i32 j = i; j <= R_ && X_[static_cast<std::size_t>(j)] != 0; ++j) {
                    const i32 len = j - i + 1;
                    if (phi + len > n) break;
                    if (m >= len) {
                        visit_E(next, phi + len, key, 't', i, j);
                    } else {
                        visit_W_interval(next, phi + m, i + m, j, j, key, 't',
                                         [&](i32, i32& x, i32& y) { x = i; y = j; });
                    }
                }
            }
        }
        // p/i-BFB: append rc(t[beta+1, gamma]) = X[n-gamma+1, n-beta].
        for (i32 gamma = 1; gamma <= phi; ++gamma) {
            const i32 start = n - gamma + 1;
            const i32 m = lce(start, tpos(phi + 1));
            const char op = gamma == phi ? 'p' : 'i';
            for (i32 beta = 0; beta < gamma; ++beta) {
                const i32 len = gamma - beta;
                if (phi + len > n) continue;
                if (m >= len) {
                    visit_E(mask, phi + len, key, op, beta, gamma);
                } else {
                    visit_W_interval(mask, phi + m, start + m, n - beta, n - beta,
                                     key, op, [&](i32, i32& x, i32& y) {
                                         x = beta;
                                         y = gamma;
                                     });
                }
            }
        }
        // f-BFB with cut c < phi: t[1,c] rc(t[beta+1, phi]).
        const i32 start = n - phi + 1;
        for (i32 c = 1; c < phi; ++c) {
            const i32 m = lce(start, tpos(c + 1));
            for (i32 beta = 0; beta < c; ++beta) {
                const i32 len = phi - beta;
                if (c + len > n) continue;
                if (m >= len) {
                    visit_E(mask, c + len, key, 'f', beta, c);
                } else {
                    visit_W_interval(mask, c + m, start + m, n - beta, n - beta,
                                     key, 'f', [&](i32, i32& x, i32& y) {
                                         x = beta;
                                         y = c;
                                     });
                }
            }
        }
    }

    void expand_W(const Record& r, u32 key) {
        const int mask = r.mask;
        const i32 n = n_;
        const i32 phi = r.phi;
        const i32 a = r.a;
        const i32 b = r.b;
        const i32 lambda = b - a + 1;
        const i32 len_s = phi + lambda;
        const i32 ra = N_ + 1 - b;  // rc(u) = X[ra, ra+lambda-1]
        const i32 tb = n - phi + 1; // rc(t[beta+1, phi]) = X[tb, n-beta]
        i32 pos = -1;

        // Case A: f-BFB with cut c <= phi; child t[1,c] rc(u) rc(t[beta+1,phi]).
        for (i32 c = 1; c <= phi; ++c) {
            const i32 M = std::min(lambda, lce(ra, tpos(c + 1)));
            if (M == lambda) {
                const i32 L2 = lce(tb, tpos(c + lambda + 1));
                // E children: beta >= phi - L2.
                const i32 beta_lo = std::max(0, phi - L2);
                visit_E_interval(mask, c + lambda + phi - (c - 1),
                                 c + lambda + phi - beta_lo, key, 'f',
                                 [&](i32 v, i32& x, i32& y) {
                                     x = c + lambda + phi - v;
                                     y = c;
                                 });
                // W children: beta <= phi - L2 - 1.
                const i32 beta_hi = std::min(c - 1, phi - L2 - 1);
                if (beta_hi >= 0) {
                    visit_W_interval(mask, c + lambda + L2, tb + L2, n - beta_hi, n,
                                     key, 'f', [&](i32 bb, i32& x, i32& y) {
                                         x = n - bb;
                                         y = c;
                                     });
                }
            } else {
                const i32 mu = lambda - M;
                const i32 v = occ(ra + M, mu, tb, pos);
                if (pos < 0) continue;
                const i32 beta_lo = std::max(0, phi - v);
                if (beta_lo > c - 1) continue;
                const i32 p = pos;
                visit_W_interval(mask, c + M, p, p + mu + phi - c,
                                 p + mu + phi - beta_lo - 1, key, 'f',
                                 [&](i32 bb, i32& x, i32& y) {
                                     x = p + mu + phi - 1 - bb;
                                     y = c;
                                 });
            }
        }

        // Case B1: f-BFB with cut c = phi + k inside the tail;
        // child t[1,phi] u[1,k] rc(s[beta+1, |s|]).
        for (i32 k = 1; k < lambda; ++k) {
            const i32 c = phi + k;
            const i32 v1 = occ(a, k, ra, pos);
            if (pos < 0) continue;
            i32 p = pos;
            i32 beta_lo;
            if (std::min(v1, lambda) < lambda) {
                beta_lo = len_s - v1;
            } else {
                const i32 v2 = occ(p, k + lambda, tb, pos);
                p = pos;
                beta_lo = std::max(0, phi - v2);
            }
            if (beta_lo > c - 1) continue;
            visit_W_interval(mask, phi, p, p + k + len_s - c,
                             p + k + len_s - beta_lo - 1, key, 'f',
                             [&](i32 bb, i32& x, i32& y) {
                                 x = p + k + len_s - 1 - bb;
                                 y = c;
                             });
        }

        // Case B2: i-BFB with q = t[beta+1, gamma], gamma <= phi;
        // child tail u rc(t[beta+1, gamma]).
        for (i32 gamma = 1; gamma <= phi; ++gamma) {
            i32 v = occ(a, lambda, n - gamma + 1, pos);
            if (pos < 0) continue;
            v = std::min(v, gamma);
            if (v == 0) continue;
            const i32 p = pos;
            visit_W_interval(mask, phi, p, p + lambda, p + lambda + v - 1, key,
                             'i', [&](i32 bb, i32& x, i32& y) {
                                 x = p + lambda + gamma - 1 - bb;
                                 y = gamma;
                             });
        }

        // Case B3: p/i-BFB with gamma = phi + g inside the tail;
        // child tail u rc(s[beta+1, gamma]).
        for (i32 g = 1; g <= lambda; ++g) {
            const i32 gamma = phi + g;
            const char op = g == lambda ? 'p' : 'i';
            const i32 s1 = N_ + 2 - a - g;  // rc(u[1,g]) = X[s1, ra+lambda-1]
            const i32 v1 = occ(a, lambda, s1, pos);
            if (pos < 0) continue;
            i32 p = pos;
            i32 beta_lo;
            if (std::min(v1, g) < g) {
                if (v1 == 0) continue;
                beta_lo = gamma - v1;
            } else {
                const i32 v2 = occ(p, lambda + g, tb, pos);
                p = pos;
                beta_lo = std::max(0, phi - v2);
            }
            visit_W_interval(mask, phi, p, p + lambda,
                             p + lambda + gamma - beta_lo - 1, key, op,
                             [&](i32 bb, i32& x, i32& y) {
                                 x = p + lambda + gamma - 1 - bb;
                                 y = gamma;
                             });
        }

        // Case C: t-BFB appending the block X[i,j]; child tail u X[i,j].
        if (may_append(mask)) {
            const int last = X_[static_cast<std::size_t>(b)];
            for (i32 i = n + 2; i <= R_; ++i) {
                if (X_[static_cast<std::size_t>(i)] == 0) continue;
                const int next = after_t(mask, last, X_[static_cast<std::size_t>(i)]);
                if (next < 0) continue;
                const i32 v = occ(a, lambda, i, pos);
                if (pos < 0 || v == 0) continue;
                const i32 p = pos;
                visit_W_interval(next, phi, p, p + lambda, p + lambda + v - 1, key,
                                 't', [&](i32 bb, i32& x, i32& y) {
                                     x = i;
                                     y = i + (bb - p - lambda);
                                 });
            }
        }
    }

    // ------------------------------------------------------------ replay

    std::vector<int> spell(const Record& r) const {
        std::vector<int> s(target_.begin(), target_.begin() + r.phi);
        if (r.kind == 'W') {
            s.insert(s.end(), X_.begin() + r.a, X_.begin() + r.b + 1);
        }
        return s;
    }

    // Rebuilds the history and checks every cycle on strings.
    std::vector<Operation> replay(u32 goal) const {
        std::vector<u32> chain;
        for (u32 id = goal; id != NO_PARENT; id = records_[id].parent) {
            chain.push_back(id);
        }
        std::reverse(chain.begin(), chain.end());
        std::vector<Operation> history;
        std::vector<int> s;
        for (std::size_t k = 1; k < chain.size(); ++k) {
            const Record& r = records_[chain[k]];
            Operation op;
            op.kind = r.op_kind;
            std::vector<int> child;
            if (r.op_kind == 't') {
                const std::vector<int> block(X_.begin() + r.op_x,
                                             X_.begin() + r.op_y + 1);
                if (r.op_x < n_ + 2 || r.op_y > R_ ||
                    std::find(block.begin(), block.end(), 0) != block.end()) {
                    throw std::logic_error("appended block is not a source");
                }
                child = s;
                child.insert(child.end(), block.begin(), block.end());
                auto found = listed_.find(block);
                op.source = found == listed_.end() ? -1 : found->second;
                op.foldback = f_of(r.mask) > f_of(records_[r.parent].mask);
            } else {
                op.x = r.op_x;
                op.y = r.op_y;
                const i32 len = static_cast<i32>(s.size());
                if (op.x < 0 || op.x >= op.y || op.y > len) {
                    throw std::logic_error("cycle out of range");
                }
                const std::vector<int> q(s.begin() + op.x, s.begin() + op.y);
                const std::vector<int> qbar = reverse_complement(q);
                if (r.op_kind == 'f') {
                    if (op.y >= len) throw std::logic_error("f-BFB without overhang");
                    child.assign(s.begin(), s.begin() + op.y);
                    const std::vector<int> rr(s.begin() + op.y, s.end());
                    const std::vector<int> rbar = reverse_complement(rr);
                    child.insert(child.end(), rbar.begin(), rbar.end());
                } else {
                    if ((r.op_kind == 'p') != (op.y == len)) {
                        throw std::logic_error("p/i-BFB kind mismatch");
                    }
                    child = s;
                }
                child.insert(child.end(), qbar.begin(), qbar.end());
            }
            if (child != spell(r)) {
                throw std::logic_error("replayed cycle does not produce the tuple");
            }
            if (r.kind == 'W') {
                const bool exact = r.phi < n_ && X_[static_cast<std::size_t>(r.a)] !=
                                                     target_[static_cast<std::size_t>(r.phi)];
                if (!exact || r.b > R_) throw std::logic_error("tuple invariant");
            }
            s = child;
            history.push_back(op);
        }
        if (s != target_) throw std::logic_error("history does not reach the target");
        return history;
    }

    // ------------------------------------------------------------ memory

    void acquire(std::size_t bytes) {
        budget_.acquire(bytes);
        owned_ += bytes;
    }
    void charge(std::size_t bytes) {
        acquire(bytes);
        search_bytes_ += bytes;
    }

    struct VecHash {
        std::size_t operator()(const std::vector<int>& v) const {
            std::size_t h = v.size();
            for (int x : v) h = h * 1000003u ^ static_cast<std::size_t>(x + 0x9e3779b9);
            return h;
        }
    };

    const std::vector<int>& target_;
    i32 n_;
    MemoryBudget& budget_;
    const JoinTable* joins_;

    std::vector<int> X_;
    i32 N_ = 0;
    i32 R_ = 0;
    int stride_ = 1;
    std::size_t width_ = 0;
    std::vector<i32> lce_;
    std::unordered_map<u64, std::pair<i32, i32>> occ_;
    std::unordered_map<std::vector<int>, int, VecHash> listed_;

    int max_t_ = 0;
    int max_fb_ = 0;
    bool count_t_ = true;
    std::vector<Record> records_;
    std::unordered_map<int, Layer> layers_;

    std::size_t owned_ = 0;
    std::size_t search_bytes_ = 0;
};

}  // namespace

Result nbfb_blocks_dist(const PreparedSources& sources,
                        const std::vector<int>& target, bool want_history,
                        MemoryBudget& budget, const JoinTable* foldback_first) {
    if (target.empty()) {
        Result result;
        result.reachable = true;
        return result;
    }
    if (sources.forward.empty()) return Result();
    BlockSearch search(sources, target, budget, foldback_first);
    return search.solve(want_history);
}

}  // namespace nbfb
