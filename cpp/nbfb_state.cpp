// Direct search over strings for the nested BFB distance; see
// nbfb_state.hpp.

#include "nbfb_state.hpp"

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

namespace nbfb {

namespace {

using Sequence = std::vector<int>;

struct SequenceHash {
    std::size_t operator()(const Sequence& sequence) const {
        std::uint64_t h = 0x9e3779b97f4a7c15ull ^ sequence.size();
        for (int value : sequence) {
            h ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(value)) +
                 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        }
        return static_cast<std::size_t>(h);
    }
};

Sequence reverse_complement(Sequence::const_iterator begin,
                            Sequence::const_iterator end) {
    Sequence result;
    result.reserve(static_cast<std::size_t>(end - begin));
    while (end != begin) {
        --end;
        result.push_back(-*end);
    }
    return result;
}

// Lexicographic cost (fold-back-shaped t-BFB joins, t-BFB cycles, total).
struct Cost {
    int foldback = 0;
    int t = 0;
    int total = 0;

    bool operator<(const Cost& other) const {
        if (foldback != other.foldback) return foldback < other.foldback;
        if (t != other.t) return t < other.t;
        return total < other.total;
    }
};

constexpr std::uint32_t NO_PARENT = 0xffffffffu;

struct Node {
    Sequence sequence;
    Cost cost;
    std::uint32_t parent = NO_PARENT;
    Operation operation;
};

// Releases the accounted bytes when the search ends, also on exceptions.
class Accounted {
  public:
    explicit Accounted(MemoryBudget& budget) : budget_(budget) {}
    ~Accounted() { budget_.release(bytes_); }
    void add(std::size_t bytes) {
        budget_.acquire(bytes);
        bytes_ += bytes;
    }

  private:
    MemoryBudget& budget_;
    std::size_t bytes_ = 0;
};

class StateSearch {
  public:
    StateSearch(const PreparedSources& sources, const Sequence& target,
                MemoryBudget& budget, const JoinTable* joins,
                bool foldable_only)
        : sources_(sources.forward),
          target_(target),
          n_(target.size()),
          joins_(joins),
          foldable_only_(foldable_only),
          accounted_(budget) {
        for (const Sequence& source : sources_) {
            if (source.empty()) {
                throw std::invalid_argument("source is empty");
            }
            source_set_.insert(source);
            source_lengths_.push_back(source.size());
        }
        std::sort(source_lengths_.begin(), source_lengths_.end());
        source_lengths_.erase(
            std::unique(source_lengths_.begin(), source_lengths_.end()),
            source_lengths_.end());

        // Every substring of the reverse complement of t: a tail u is
        // foldable exactly when u occurs in it.
        const Sequence tbar = reverse_complement(target_.begin(), target_.end());
        accounted_.add(2 * n_ * n_ * (sizeof(int) * (n_ + 2) / 2 + 48));
        for (std::size_t a = 0; a < n_; ++a) {
            for (std::size_t b = a + 1; b <= n_; ++b) {
                window_set_.emplace(tbar.begin() + a, tbar.begin() + b);
            }
        }
        // With foldable_only, a tail may also occur in t itself.
        if (foldable_only_) {
            accounted_.add(n_ * n_ * (sizeof(int) * (n_ + 2) / 2 + 48));
            for (std::size_t a = 0; a < n_; ++a) {
                for (std::size_t b = a + 1; b <= n_; ++b) {
                    window_set_.emplace(target_.begin() + a, target_.begin() + b);
                }
            }
        }
    }

    Result solve(bool want_history) {
        std::vector<std::vector<std::uint32_t>> by_length(n_ + 1);
        add_node(Sequence(), Cost(), NO_PARENT, Operation(), by_length);

        for (std::size_t length = 0; length < n_; ++length) {
            // Nodes of this length are final: every edge into them comes
            // from a strictly shorter string, which has been expanded.
            for (std::size_t k = 0; k < by_length[length].size(); ++k) {
                expand(by_length[length][k], by_length);
            }
        }

        Result result;
        auto found = index_.find(target_);
        if (found == index_.end()) return result;
        const Node& end = nodes_[found->second];
        result.reachable = true;
        result.foldback_t = end.cost.foldback;
        result.t_cycles = end.cost.t;
        result.total = end.cost.total;
        if (want_history) {
            for (std::uint32_t id = found->second; nodes_[id].parent != NO_PARENT;
                 id = nodes_[id].parent) {
                result.history.push_back(nodes_[id].operation);
            }
            std::reverse(result.history.begin(), result.history.end());
        }
        return result;
    }

  private:
    std::size_t lcp(const Sequence& sequence) const {
        std::size_t k = 0;
        while (k < sequence.size() && k < n_ && sequence[k] == target_[k]) ++k;
        return k;
    }

    // Tail beyond the longest common prefix is empty or the reverse
    // complement of a substring of t.
    bool foldable(const Sequence& sequence, std::size_t phi) const {
        if (phi == sequence.size()) return true;
        return window_set_.count(Sequence(sequence.begin() + phi, sequence.end())) > 0;
    }

    // t[1,psi] y with y empty, an inverted window of t, or a source. A
    // suffix of an inverted window is an inverted window, so the first two
    // cases are tested at psi = phi.
    bool nice(const Sequence& sequence, std::size_t phi) const {
        if (foldable(sequence, phi)) return true;
        if (foldable_only_) return false;
        for (std::size_t length : source_lengths_) {
            if (length > sequence.size()) break;
            if (sequence.size() - length > phi) continue;
            if (source_set_.count(Sequence(sequence.end() - length, sequence.end()))) {
                return true;
            }
        }
        return false;
    }

    bool admitted(const Sequence& child) const { return nice(child, lcp(child)); }

    void add_node(Sequence sequence, const Cost& cost, std::uint32_t parent,
                  const Operation& operation,
                  std::vector<std::vector<std::uint32_t>>& by_length) {
        auto found = index_.find(sequence);
        if (found != index_.end()) {
            Node& node = nodes_[found->second];
            if (cost < node.cost) {
                node.cost = cost;
                node.parent = parent;
                node.operation = operation;
            }
            return;
        }
        accounted_.add(2 * sizeof(int) * sequence.size() + 160);
        const std::uint32_t id = static_cast<std::uint32_t>(nodes_.size());
        by_length[sequence.size()].push_back(id);
        index_.emplace(sequence, id);
        Node node;
        node.sequence = std::move(sequence);
        node.cost = cost;
        node.parent = parent;
        node.operation = operation;
        nodes_.push_back(std::move(node));
    }

    void expand(std::uint32_t id,
                std::vector<std::vector<std::uint32_t>>& by_length) {
        // Copies: add_node may reallocate nodes_.
        const Sequence s = nodes_[id].sequence;
        const Cost cost = nodes_[id].cost;
        const std::size_t length = s.size();

        Cost fold = cost;
        fold.total += 1;
        for (std::size_t x = 0; x < length; ++x) {
            for (std::size_t y = x + 1; y <= length; ++y) {
                if (length + (y - x) > n_) break;
                const Sequence qbar = reverse_complement(s.begin() + x, s.begin() + y);
                Operation operation;
                operation.x = static_cast<int>(x);
                operation.y = static_cast<int>(y);

                // p (R empty) or i: P Q R -> P Q R Qbar.
                Sequence keep(s);
                keep.insert(keep.end(), qbar.begin(), qbar.end());
                operation.kind = y == length ? 'p' : 'i';
                if (admitted(keep)) {
                    add_node(std::move(keep), fold, id, operation, by_length);
                }

                // f: P Q R -> P Q Rbar Qbar.
                if (y < length) {
                    Sequence flip(s.begin(), s.begin() + y);
                    const Sequence rbar = reverse_complement(s.begin() + y, s.end());
                    flip.insert(flip.end(), rbar.begin(), rbar.end());
                    flip.insert(flip.end(), qbar.begin(), qbar.end());
                    operation.kind = 'f';
                    if (admitted(flip)) {
                        add_node(std::move(flip), fold, id, operation, by_length);
                    }
                }
            }
        }

        for (std::size_t k = 0; k < sources_.size(); ++k) {
            const Sequence& source = sources_[k];
            if (length + source.size() > n_) continue;
            Sequence child(s);
            child.insert(child.end(), source.begin(), source.end());
            if (!admitted(child)) continue;
            Operation operation;
            operation.kind = 't';
            operation.source = static_cast<int>(k);
            operation.foldback = joins_ != nullptr && length > 0 &&
                                 joins_->is_foldback(s.back(), source.front());
            Cost append = cost;
            append.t += 1;
            append.total += 1;
            if (operation.foldback) append.foldback += 1;
            add_node(std::move(child), append, id, operation, by_length);
        }
    }

    const std::vector<Sequence>& sources_;
    const Sequence& target_;
    std::size_t n_;
    const JoinTable* joins_;
    bool foldable_only_;
    Accounted accounted_;

    std::unordered_set<Sequence, SequenceHash> source_set_;
    std::vector<std::size_t> source_lengths_;
    std::unordered_set<Sequence, SequenceHash> window_set_;

    std::vector<Node> nodes_;
    std::unordered_map<Sequence, std::uint32_t, SequenceHash> index_;
};

}  // namespace

Result nbfb_state_dist(const PreparedSources& sources,
                       const std::vector<int>& target,
                       bool want_history,
                       MemoryBudget& budget,
                       const JoinTable* foldback_first,
                       bool foldable_only) {
    if (target.empty()) {
        Result result;
        result.reachable = true;
        return result;
    }
    StateSearch search(sources, target, budget, foldback_first, foldable_only);
    return search.solve(want_history);
}

}  // namespace nbfb
