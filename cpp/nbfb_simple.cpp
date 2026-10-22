// Nested BFB distance over a general source set; see nbfb_simple.hpp.

#include "nbfb_simple.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

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

// ------------------------------------------------------------ hashing

// Polynomial hashing modulo the Mersenne prime 2^61 - 1. Every lookup is
// confirmed by a direct comparison, so a collision never changes a result.
constexpr std::uint64_t MOD = (std::uint64_t(1) << 61) - 1;

std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
    const unsigned __int128 product = static_cast<unsigned __int128>(a) * b;
    const std::uint64_t low = static_cast<std::uint64_t>(product & MOD);
    const std::uint64_t high = static_cast<std::uint64_t>(product >> 61);
    std::uint64_t sum = low + high;
    if (sum >= MOD) sum -= MOD;
    return sum;
}

std::uint64_t add(std::uint64_t a, std::uint64_t b) {
    std::uint64_t sum = a + b;
    if (sum >= MOD) sum -= MOD;
    return sum;
}

std::uint64_t sub(std::uint64_t a, std::uint64_t b) {
    return a >= b ? a - b : a + MOD - b;
}

std::uint64_t symbol_value(int symbol) {
    return static_cast<std::uint64_t>(static_cast<std::int64_t>(symbol) +
                                      (std::int64_t(1) << 32));
}

// Key of a string of a given length with a given hash.
std::uint64_t window_key(std::uint64_t hash, std::uint32_t length) {
    return hash ^ (static_cast<std::uint64_t>(length) * 0x9e3779b97f4a7c15ull);
}

// ------------------------------------------------------------ costs

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
    bool operator==(const Cost& other) const {
        return foldback == other.foldback && t == other.t &&
               total == other.total;
    }
};

const Cost INFINITE_COST = {std::numeric_limits<int>::max(), 0, 0};

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

// Open-addressing hash map from 64-bit keys to 32-bit values, with linear
// probing. Its tables are accounted in the search's budget.
class FlatMap {
  public:
    static constexpr std::uint32_t NONE = 0xffffffffu;

    explicit FlatMap(Accounted& accounted) : accounted_(accounted) {}

    std::uint32_t find(std::uint64_t key) const {
        if (values_.empty()) return NONE;
        for (std::size_t slot = mix(key) & mask_;; slot = (slot + 1) & mask_) {
            if (values_[slot] == NONE) return NONE;
            if (keys_[slot] == key) return values_[slot];
        }
    }

    // The key must be absent.
    void insert(std::uint64_t key, std::uint32_t value) {
        if (2 * (size_ + 1) > values_.size()) grow();
        std::size_t slot = mix(key) & mask_;
        while (values_[slot] != NONE) slot = (slot + 1) & mask_;
        keys_[slot] = key;
        values_[slot] = value;
        ++size_;
    }

    void clear() {
        std::fill(values_.begin(), values_.end(), NONE);
        size_ = 0;
    }

  private:
    static std::size_t mix(std::uint64_t key) {
        key ^= key >> 33;
        key *= 0xff51afd7ed558ccdull;
        key ^= key >> 33;
        return static_cast<std::size_t>(key);
    }

    void grow() {
        const std::size_t capacity = values_.empty() ? 1024 : 2 * values_.size();
        accounted_.add((capacity - values_.size()) *
                       (sizeof(std::uint64_t) + sizeof(std::uint32_t)));
        std::vector<std::uint64_t> keys(capacity);
        std::vector<std::uint32_t> values(capacity, NONE);
        keys_.swap(keys);
        values_.swap(values);
        mask_ = capacity - 1;
        size_ = 0;
        for (std::size_t slot = 0; slot < values.size(); ++slot) {
            if (values[slot] != NONE) insert(keys[slot], values[slot]);
        }
    }

    Accounted& accounted_;
    std::vector<std::uint64_t> keys_;
    std::vector<std::uint32_t> values_;
    std::size_t mask_ = 0;
    std::size_t size_ = 0;
};

// ------------------------------------------------------------ ropes

// A piece T[off, off+len) of the text T = t # rc(t).
struct Piece {
    std::uint32_t off;
    std::uint32_t len;
};

// A string given as a concatenation of at most 8 pieces of T.
struct Rope {
    Piece piece[8];
    int count = 0;
    std::uint32_t length = 0;

    void push(std::uint32_t off, std::uint32_t len) {
        if (len == 0) return;
        piece[count++] = {off, len};
        length += len;
    }
    void append(const Rope& other) {
        for (int k = 0; k < other.count; ++k) {
            push(other.piece[k].off, other.piece[k].len);
        }
    }
};

constexpr std::uint32_t NO_PARENT = 0xffffffffu;
constexpr std::int32_t NOT_FOLDED = -2;

// Cost classes of chain edges: all edges of one chain add the same cost,
// so the first offer of a class to a node is the best offer of the class.
constexpr int FOLD = 0;           // p-, i- and f-cycles: (0, 0, 1)
constexpr int APPEND = 1;         // t-cycles: (0, 1, 1)
constexpr int APPEND_FOLDBACK = 2;  // fold-back-shaped t-cycles: (1, 1, 1)

struct Node {
    std::uint32_t phi;   // longest common prefix with t
    std::uint32_t off;   // first occurrence of the tail in T
    std::uint32_t len;   // length of the tail
    Cost cost = INFINITE_COST;
    std::uint32_t parent = NO_PARENT;
    Operation operation;
    // Per cost class: NOT_FOLDED until a chain edge of the class reaches
    // the node. Afterwards every prefix of the node whose length lies
    // strictly between jump and the length of the node has been reached by
    // a chain edge of the class.
    std::int32_t jump[3] = {NOT_FOLDED, NOT_FOLDED, NOT_FOLDED};
    bool settled = false;
};

class SimpleSearch {
  public:
    SimpleSearch(const PreparedSources& sources, const Sequence& target,
                 MemoryBudget& budget, const JoinTable* joins, bool blocks)
        : sources_(sources.forward),
          target_(target),
          n_(static_cast<std::uint32_t>(target.size())),
          joins_(joins),
          blocks_(blocks),
          accounted_(budget),
          windows_(accounted_),
          root_(accounted_),
          bigrams_(accounted_),
          follow_index_(accounted_),
          index_(accounted_) {
        if (target.size() >= (std::size_t(1) << 20)) {
            throw std::invalid_argument("target too long");
        }
        for (const Sequence& source : sources_) {
            if (source.empty()) throw std::invalid_argument("source is empty");
        }
        build_text();
        build_lce();
        build_windows();
        if (blocks_) {
            build_blocks();
        } else {
            build_sources();
        }
        build_extensions();
    }

    Result solve(bool want_history) {
        const std::uint32_t start = node_id(0, 0, 0);
        const std::uint32_t goal = node_id(n_, 0, 0);
        offer(start, Cost{0, 0, 0}, NO_PARENT, Operation());

        while (!heap_.empty()) {
            std::pop_heap(heap_.begin(), heap_.end(), Later());
            const auto [cost, id] = heap_.back();
            heap_.pop_back();
            if (nodes_[id].settled || !(nodes_[id].cost == cost)) continue;
            nodes_[id].settled = true;
            if (id == goal) break;
            expand(id);
        }

        Result result;
        const Node& end = nodes_[goal];
        if (!end.settled) return result;
        result.reachable = true;
        result.foldback_t = end.cost.foldback;
        result.t_cycles = end.cost.t;
        result.total = end.cost.total;
        for (std::uint32_t id = goal; nodes_[id].parent != NO_PARENT;
             id = nodes_[id].parent) {
            result.history.push_back(nodes_[id].operation);
        }
        std::reverse(result.history.begin(), result.history.end());
        replay(result);
        if (blocks_) {
            // Index of the appended block among the source lines, or -1.
            for (Operation& op : result.history) {
                if (op.kind != 't') continue;
                const Sequence block(text_.begin() + op.x, text_.begin() + op.x + op.y);
                auto found = listed_.find(block);
                op.source = found == listed_.end() ? -1 : found->second;
                op.x = -1;
                op.y = -1;
            }
        }
        if (!want_history) result.history.clear();
        return result;
    }

  private:
    // -------------------------------------------------------- tables

    void build_text() {
        // T = t # rc(t), with # = 0 (segment identifiers are nonzero).
        text_.reserve(2 * n_ + 1);
        text_.insert(text_.end(), target_.begin(), target_.end());
        text_.push_back(0);
        for (std::uint32_t a = n_; a > 0; --a) text_.push_back(-target_[a - 1]);
        if (blocks_) {
            // T = t # rc(t) # g_1 # rc(g_1) # ... # g_r # rc(g_r) for the
            // maximal source lines g_1, ..., g_r.
            for (const Sequence& g : maximal_sources()) {
                const Sequence g_bar = reverse_complement(g);
                for (const Sequence* line : {&g, &g_bar}) {
                    text_.push_back(0);
                    const std::uint32_t begin = static_cast<std::uint32_t>(text_.size());
                    text_.insert(text_.end(), line->begin(), line->end());
                    strings_.push_back({begin, static_cast<std::uint32_t>(text_.size())});
                }
            }
        }
        accounted_.add(text_.size() * (sizeof(int) + 2 * sizeof(std::uint64_t)));
    }

    static Sequence reverse_complement(const Sequence& s) {
        Sequence r(s.rbegin(), s.rend());
        for (int& v : r) v = -v;
        return r;
    }

    static bool occurs(const Sequence& text, const Sequence& pattern) {
        return std::search(text.begin(), text.end(), pattern.begin(), pattern.end()) !=
               text.end();
    }

    // The source lines, longest first, without every line that occurs in a
    // kept line or in the reverse complement of a kept line. Their blocks
    // and the reverse complements of the blocks are the source set.
    std::vector<Sequence> maximal_sources() const {
        std::vector<Sequence> all = sources_;
        std::stable_sort(all.begin(), all.end(), [](const Sequence& l, const Sequence& r) {
            return l.size() > r.size();
        });
        std::vector<Sequence> kept;
        for (const Sequence& x : all) {
            bool inside = false;
            for (const Sequence& g : kept) {
                if (occurs(g, x) || occurs(reverse_complement(g), x)) {
                    inside = true;
                    break;
                }
            }
            if (!inside) kept.push_back(x);
        }
        return kept;
    }

    void build_blocks() {
        std::size_t bytes = 0;
        for (std::size_t k = 0; k < sources_.size(); ++k) {
            bytes += sources_[k].size() * sizeof(int) + 64;
            listed_.emplace(sources_[k], static_cast<int>(k));
        }
        accounted_.add(bytes);
    }

    void build_lce() {
        // lce_[i * width + j] = longest common extension of T[i..] and T[j..].
        width_ = text_.size() + 1;
        if (width_ > 65535) throw std::invalid_argument("target too long");
        accounted_.add(width_ * width_ * sizeof(std::uint16_t));
        lce_.assign(width_ * width_, 0);
        for (std::size_t i = text_.size(); i-- > 0;) {
            for (std::size_t j = text_.size(); j-- > 0;) {
                if (text_[i] == text_[j] && text_[i] != 0) {
                    lce_[i * width_ + j] = static_cast<std::uint16_t>(
                        lce_[(i + 1) * width_ + (j + 1)] + 1);
                }
            }
        }
    }

    std::uint32_t lce(std::uint32_t i, std::uint32_t j) const {
        return lce_[i * width_ + j];
    }

    void build_windows() {
        // Every window of t and of rc(t), keyed by hash and length, mapped
        // to its first occurrence in T. A collision between two different
        // windows restarts with another base.
        for (std::uint64_t seed = 1;; ++seed) {
            base_ = add(mul(0x1234567890abcdefull % MOD, seed),
                        1000003 + 7919 * seed) % MOD;
            if (base_ < 256) continue;
            hash_text();
            windows_.clear();
            bool collision = false;
            for (std::uint32_t len = 1; len <= n_ && !collision; ++len) {
                for (std::uint32_t off = 0; off < text_.size(); ++off) {
                    if (!inside(off, len)) continue;
                    const std::uint64_t key = window_key(piece_hash(off, len), len);
                    const std::uint32_t found = windows_.find(key);
                    if (found == FlatMap::NONE) {
                        windows_.insert(key, off);
                    } else if (lce(found, off) < len) {
                        collision = true;
                        break;
                    }
                }
            }
            if (!collision) break;
        }
    }

    // T[off, off+len) lies inside t or inside rc(t).
    bool inside(std::uint32_t off, std::uint32_t len) const {
        if (off < n_) return off + len <= n_;
        return off > n_ && off + len <= 2 * n_ + 1;
    }

    void hash_text() {
        prefix_.assign(text_.size() + 1, 0);
        const std::size_t longest = text_.size();
        power_.assign(longest + 1, 1);
        for (std::size_t a = 1; a <= longest; ++a) power_[a] = mul(power_[a - 1], base_);
        for (std::size_t a = 0; a < text_.size(); ++a) {
            prefix_[a + 1] = add(mul(prefix_[a], base_), symbol_value(text_[a]));
        }
    }

    std::uint64_t piece_hash(std::uint32_t off, std::uint32_t len) const {
        return sub(prefix_[off + len], mul(prefix_[off], power_[len]));
    }

    void build_sources() {
        // A trie of the sources. Every trie node is a non-empty prefix of a
        // source, and lists the sources that end there.
        std::size_t total = 0;
        for (const Sequence& source : sources_) total += source.size();
        accounted_.add((total + 1) * (sizeof(TrieNode) + 16));
        trie_.emplace_back();
        for (std::size_t k = 0; k < sources_.size(); ++k) {
            std::uint32_t at = 0;
            for (int symbol : sources_[k]) {
                std::uint32_t next = NO_PARENT;
                for (const auto& [label, child] : trie_[at].children) {
                    if (label == symbol) next = child;
                }
                if (next == NO_PARENT) {
                    next = static_cast<std::uint32_t>(trie_.size());
                    trie_[at].children.emplace_back(symbol, next);
                    trie_.emplace_back();
                }
                at = next;
            }
            trie_[at].sources.push_back(static_cast<std::uint32_t>(k));
        }
        for (const auto& [symbol, child] : trie_[0].children) {
            root_.insert(symbol_key(symbol), child);
        }
    }

    static std::uint64_t symbol_key(int symbol) {
        return static_cast<std::uint64_t>(static_cast<std::uint32_t>(symbol));
    }

    static std::uint64_t pair_key(int left, int right) {
        return (symbol_key(left) << 32) | symbol_key(right);
    }

    void build_extensions() {
        // The adjacent pairs of t and of rc(t), and for every window of T
        // (by first occurrence) the symbols that follow some occurrence.
        int largest = 0;
        for (int symbol : text_) largest = std::max(largest, std::abs(symbol));
        symbols_ = largest;
        accounted_.add((2 * static_cast<std::size_t>(largest) + 1) * sizeof(std::uint32_t));
        stamp_.assign(2 * static_cast<std::size_t>(largest) + 1, 0);

        std::vector<std::pair<std::uint64_t, int>> follow;
        accounted_.add(static_cast<std::size_t>(n_) * n_ *
                       (sizeof(std::pair<std::uint64_t, int>) + sizeof(int)));
        for (std::uint32_t off = 0; off + 1 < text_.size(); ++off) {
            if (!inside(off, 2)) continue;
            if (bigrams_.find(pair_key(text_[off], text_[off + 1])) == FlatMap::NONE) {
                bigrams_.insert(pair_key(text_[off], text_[off + 1]), 1);
            }
            for (std::uint32_t len = 1; inside(off, len + 1); ++len) {
                const std::uint32_t first =
                    windows_.find(window_key(piece_hash(off, len), len));
                follow.emplace_back(key_of(0, first, len), text_[off + len]);
            }
        }
        std::sort(follow.begin(), follow.end());
        follow.erase(std::unique(follow.begin(), follow.end()), follow.end());
        for (std::size_t k = 0; k < follow.size(); ++k) {
            if (k == 0 || follow[k].first != follow[k - 1].first) {
                follow_index_.insert(follow[k].first,
                                     static_cast<std::uint32_t>(follow_begin_.size()));
                follow_begin_.push_back(static_cast<std::uint32_t>(k));
            }
            follow_symbols_.push_back(follow[k].second);
        }
        follow_begin_.push_back(static_cast<std::uint32_t>(follow.size()));
    }

    // -------------------------------------------------------- ropes

    // Reverse complement of a piece: the mirrored piece of T.
    Piece mirror(Piece p) const { return {2 * n_ + 1 - p.off - p.len, p.len}; }

    Rope slice(const Rope& rope, std::uint32_t from, std::uint32_t to) const {
        Rope result;
        std::uint32_t at = 0;
        for (int k = 0; k < rope.count && at < to; ++k) {
            const Piece& p = rope.piece[k];
            const std::uint32_t begin = std::max(from, at);
            const std::uint32_t end = std::min(to, at + p.len);
            if (begin < end) result.push(p.off + (begin - at), end - begin);
            at += p.len;
        }
        return result;
    }

    Rope reversed(const Rope& rope) const {
        Rope result;
        for (int k = rope.count; k-- > 0;) {
            const Piece p = mirror(rope.piece[k]);
            result.push(p.off, p.len);
        }
        return result;
    }

    // Longest common prefix of a rope and t.
    std::uint32_t prefix_match(const Rope& rope) const {
        std::uint32_t p = 0;
        for (int k = 0; k < rope.count && p < n_; ++k) {
            const Piece& piece = rope.piece[k];
            const std::uint32_t m =
                std::min({piece.len, n_ - p, lce(piece.off, p)});
            p += m;
            if (m < piece.len) break;
        }
        return p;
    }

    std::uint64_t rope_hash(const Rope& rope) const {
        std::uint64_t h = 0;
        for (int k = 0; k < rope.count; ++k) {
            const Piece& p = rope.piece[k];
            h = add(mul(h, power_[p.len]), piece_hash(p.off, p.len));
        }
        return h;
    }

    // First occurrence in T of a non-empty rope, or -1.
    std::int64_t find_window(const Rope& rope) const {
        if (rope.length > n_) return -1;
        const std::uint32_t found =
            windows_.find(window_key(rope_hash(rope), rope.length));
        if (found == FlatMap::NONE) return -1;
        std::uint32_t at = found;
        for (int k = 0; k < rope.count; ++k) {
            const Piece& p = rope.piece[k];
            if (lce(p.off, at) < p.len) return -1;
            at += p.len;
        }
        return found;
    }

    // -------------------------------------------------------- nodes

    static std::uint64_t key_of(std::uint32_t phi, std::uint32_t off,
                                std::uint32_t len) {
        return (static_cast<std::uint64_t>(phi) << 42) |
               (static_cast<std::uint64_t>(off) << 21) | len;
    }

    std::uint32_t node_id(std::uint32_t phi, std::uint32_t off,
                          std::uint32_t len) {
        if (len == 0) off = 0;
        const std::uint64_t key = key_of(phi, off, len);
        const std::uint32_t found = index_.find(key);
        if (found != FlatMap::NONE) return found;
        reserve(nodes_);
        const std::uint32_t id = static_cast<std::uint32_t>(nodes_.size());
        Node node;
        node.phi = phi;
        node.off = off;
        node.len = len;
        nodes_.push_back(node);
        index_.insert(key, id);
        return id;
    }

    void offer(std::uint32_t id, const Cost& cost, std::uint32_t parent,
               const Operation& operation) {
        Node& node = nodes_[id];
        if (node.settled || !(cost < node.cost)) return;
        node.cost = cost;
        node.parent = parent;
        node.operation = operation;
        reserve(heap_);
        heap_.push_back({cost, id});
        std::push_heap(heap_.begin(), heap_.end(), Later());
    }

    // Grows a vector by doubling before a push_back, accounting the new
    // capacity in the budget.
    template <typename T>
    void reserve(std::vector<T>& items) {
        if (items.size() < items.capacity()) return;
        const std::size_t capacity = std::max<std::size_t>(1024, 2 * items.capacity());
        accounted_.add((capacity - items.capacity()) * sizeof(T));
        items.reserve(capacity);
    }

    // The node of the prefix X[0, length) of a rope X with longest common
    // prefix phi_x with t, or -1 if that prefix is not a foldable prefix.
    std::int64_t prefix_node(const Rope& x, std::uint32_t phi_x,
                             std::uint32_t length, bool create) {
        const std::uint32_t phi = std::min(phi_x, length);
        if (phi == length) return create ? node_id(phi, 0, 0) : 0;
        const std::int64_t off = find_window(slice(x, phi, length));
        if (off < 0) return -1;
        if (!create) return 0;
        return node_id(phi, static_cast<std::uint32_t>(off), length - phi);
    }

    Rope rope_of(const Node& node) const {
        Rope rope;
        rope.push(0, node.phi);
        rope.push(node.off, node.len);
        return rope;
    }

    int last_symbol(const Node& node) const {
        if (node.len > 0) return text_[node.off + node.len - 1];
        return target_[node.phi - 1];
    }

    // A foldable prefix t[0, phi) T[off, off+len) and the hash of its tail,
    // extended one symbol at a time.
    struct Extension {
        std::uint32_t phi;
        std::uint32_t off;
        std::uint32_t len;
        std::uint64_t hash;
    };

    Extension extension_of(const Node& node) const {
        return {node.phi, node.off, node.len, piece_hash(node.off, node.len)};
    }

    // Appends a symbol; false if the result is not a foldable prefix.
    bool extend(Extension& e, int symbol) const {
        if (e.len == 0 && e.phi < n_ && target_[e.phi] == symbol) {
            ++e.phi;
            return true;
        }
        if (e.phi + e.len + 1 > n_) return false;
        const std::uint64_t hash = add(mul(e.hash, base_), symbol_value(symbol));
        const std::uint32_t off = windows_.find(window_key(hash, e.len + 1));
        if (off == FlatMap::NONE) return false;
        if (e.len > 0 && lce(e.off, off) < e.len) return false;
        if (text_[off + e.len] != symbol) return false;
        e.off = off;
        e.len += 1;
        e.hash = hash;
        return true;
    }

    // -------------------------------------------------------- search

    std::int32_t symbol_at(const Node& node, std::uint32_t a) const {
        return a < node.phi ? target_[a] : text_[node.off + (a - node.phi)];
    }

    void expand(std::uint32_t id) {
        const Node u = nodes_[id];  // copy: nodes_ may reallocate
        const std::uint32_t length = u.phi + u.len;
        const Rope w = rope_of(u);
        const Extension base = extension_of(u);

        // The symbols c with w c a foldable prefix: every symbol if w is a
        // prefix of t, else the symbols that follow the tail in T.
        std::uint32_t follow_from = 0;
        std::uint32_t follow_to = 0;
        if (u.len > 0) {
            const std::uint32_t k = follow_index_.find(key_of(0, u.off, u.len));
            if (k != FlatMap::NONE) {
                follow_from = follow_begin_[k];
                follow_to = follow_begin_[k + 1];
            }
            ++epoch_;
            for (std::uint32_t k2 = follow_from; k2 < follow_to; ++k2) {
                stamp_[follow_symbols_[k2] + symbols_] = epoch_;
            }
        }
        auto follows = [&](int symbol) {
            return u.len == 0 || stamp_[symbol + symbols_] == epoch_;
        };

        // The shortest suffix start j0 with w[j0, length) a window of T.
        std::uint32_t j0 = 0;
        if (u.len > 0 && length < n_) {
            std::uint32_t good = u.phi;  // w[phi, length) is the tail
            std::uint32_t bad = 0;
            if (find_window(slice(w, 0, length)) >= 0) {
                good = 0;
            } else {
                while (good - bad > 1) {
                    const std::uint32_t mid = bad + (good - bad) / 2;
                    if (find_window(slice(w, mid, length)) >= 0) good = mid;
                    else bad = mid;
                }
            }
            j0 = good;
        }

        Cost fold = u.cost;
        fold.total += 1;
        for (std::uint32_t b = 1; b <= length && length < n_; ++b) {
            const std::uint32_t high = std::min(length + b, n_);
            const int before = symbol_at(u, b - 1);
            Rope head;
            Rope head_bar;
            bool sliced = false;
            auto slice_head = [&]() {
                if (sliced) return;
                head = slice(w, 0, b);
                head_bar = reversed(head);
                sliced = true;
            };

            // p (b = length) or i: P Q R -> P Q R Qbar with Q = w[a, b).
            // The shortest child is w followed by -w[b-1].
            if (follows(-before)) {
                slice_head();
                Rope keep = w;
                keep.append(head_bar);
                const char kind = b == length ? 'p' : 'i';
                chain(keep, length, high, FOLD, id, fold,
                      [&](std::uint32_t l) { return fold_operation(kind, b, l - length); });
            }

            // f: P Q R -> P Q Rbar Qbar with Q = w[a, b), R = w[b, length).
            // The shortest child is w[0, b) rc(w[b-1, length)). If it is a
            // foldable prefix, then the pair at positions b-1, b is
            // adjacent in T, or the child's common prefix with t is
            // w[0, b) and rc(w[b-1, length)) is a window of T.
            if (b < length) {
                const int last = symbol_at(u, length - 1);
                const bool adjacent =
                    bigrams_.find(pair_key(before, -last)) != FlatMap::NONE;
                if (adjacent || (b <= u.phi && b - 1 >= j0)) {
                    slice_head();
                    Rope flip = head;
                    flip.append(reversed(slice(w, b, length)));
                    flip.append(head_bar);
                    chain(flip, length, high, FOLD, id, fold,
                          [&](std::uint32_t l) { return fold_operation('f', b, l - length); });
                }
            }
        }

        // t: walk the trie of the sources, pruning at the first prefix that
        // is not a foldable prefix (foldable prefixes are closed under
        // taking prefixes).
        auto start = [&](int symbol, std::uint32_t child) {
            Extension e = base;
            if (!extend(e, symbol)) return;
            Operation operation;
            operation.kind = 't';
            operation.foldback = joins_ != nullptr && length > 0 &&
                                 joins_->is_foldback(last_symbol(u), symbol);
            Cost cost = u.cost;
            cost.t += 1;
            cost.total += 1;
            if (operation.foldback) cost.foldback += 1;
            walk(child, e, id, operation, cost);
        };
        if (blocks_) {
            // Every block starting at position c of a string g of T is a
            // prefix of the rest of g, so the children s x for these blocks
            // are the prefixes of s g[c..] longer than s: one chain per
            // start position.
            const int last = length > 0 ? last_symbol(u) : 0;
            for (const auto& [begin, end] : strings_) {
                for (std::uint32_t c = begin; c < end && length < n_; ++c) {
                    const int symbol = text_[c];
                    if (std::abs(symbol) > symbols_ || !follows(symbol)) continue;
                    const bool foldback = joins_ != nullptr && length > 0 &&
                                          joins_->is_foldback(last, symbol);
                    Cost cost = u.cost;
                    cost.t += 1;
                    cost.total += 1;
                    if (foldback) cost.foldback += 1;
                    Rope x = w;
                    x.push(c, end - c);
                    const std::uint32_t high = std::min(length + (end - c), n_);
                    chain(x, length, high, foldback ? APPEND_FOLDBACK : APPEND, id, cost,
                          [&](std::uint32_t l) {
                              Operation operation;
                              operation.kind = 't';
                              operation.foldback = foldback;
                              operation.x = static_cast<int>(c);
                              operation.y = static_cast<int>(l - length);
                              return operation;
                          });
                }
            }
        } else if (u.len == 0) {
            for (const auto& [symbol, child] : trie_[0].children) start(symbol, child);
        } else {
            for (std::uint32_t k = follow_from; k < follow_to; ++k) {
                const int symbol = follow_symbols_[k];
                const std::uint32_t child = root_.find(symbol_key(symbol));
                if (child != FlatMap::NONE) start(symbol, child);
            }
        }
    }

    void walk(std::uint32_t at, const Extension& e, std::uint32_t parent,
              Operation& operation, const Cost& cost) {
        for (std::uint32_t k : trie_[at].sources) {
            operation.source = static_cast<int>(k);
            offer(node_id(e.phi, e.off, e.len), cost, parent, operation);
        }
        for (const auto& [symbol, child] : trie_[at].children) {
            Extension next = e;
            if (extend(next, symbol)) walk(child, next, parent, operation, cost);
        }
    }

    // The cycle P Q R -> ... with |P Q| = b and |Q| = gain.
    static Operation fold_operation(char kind, std::uint32_t b, std::uint32_t gain) {
        Operation operation;
        operation.kind = kind;
        operation.x = static_cast<int>(b - gain);
        operation.y = static_cast<int>(b);
        return operation;
    }

    // Offers the children x[0, l) for l in (length, high] of node parent,
    // all with the given cost of the given class; make(l) is the cycle that
    // produces x[0, l).
    template <typename Make>
    void chain(const Rope& x, std::uint32_t length, std::uint32_t high, int cls,
               std::uint32_t parent, const Cost& cost, const Make& make) {
        const std::uint32_t phi_x = prefix_match(slice(x, 0, high));
        const std::uint32_t low = length + 1;
        if (prefix_node(x, phi_x, low, false) < 0) return;
        // Largest foldable length; foldable prefixes are closed under
        // taking prefixes.
        std::uint32_t good = low;
        std::uint32_t bad = high + 1;
        if (prefix_node(x, phi_x, high, false) >= 0) {
            good = high;
        } else {
            bad = high;
            while (bad - good > 1) {
                const std::uint32_t mid = good + (bad - good) / 2;
                if (prefix_node(x, phi_x, mid, false) >= 0) good = mid;
                else bad = mid;
            }
        }

        path_.clear();
        std::int64_t at = good;
        while (at >= static_cast<std::int64_t>(low)) {
            const std::uint32_t l = static_cast<std::uint32_t>(at);
            const std::uint32_t child =
                static_cast<std::uint32_t>(prefix_node(x, phi_x, l, true));
            path_.push_back(child);
            if (nodes_[child].jump[cls] == NOT_FOLDED) {
                nodes_[child].jump[cls] = static_cast<std::int32_t>(l) - 1;
                offer(child, cost, parent, make(l));
                at = l - 1;
            } else {
                at = nodes_[child].jump[cls];
            }
        }
        for (std::uint32_t child : path_) {
            nodes_[child].jump[cls] = static_cast<std::int32_t>(at);
        }
    }

    // -------------------------------------------------------- checking

    // Replays the history on strings: every cycle must be valid, every
    // string a foldable prefix of t, the last string t, and the cost the
    // reported one.
    void replay(const Result& result) const {
        Sequence s;
        Cost cost;
        for (const Operation& op : result.history) {
            const std::size_t length = s.size();
            if (op.kind == 't') {
                const Sequence x = appended(op);
                const bool foldback = joins_ != nullptr && length > 0 &&
                                      joins_->is_foldback(s.back(), x.front());
                if (foldback != op.foldback) fail("fold-back flag");
                s.insert(s.end(), x.begin(), x.end());
                cost.t += 1;
                if (foldback) cost.foldback += 1;
            } else {
                const std::size_t a = static_cast<std::size_t>(op.x);
                const std::size_t b = static_cast<std::size_t>(op.y);
                if (!(a < b && b <= length)) fail("cut points");
                if ((op.kind == 'p') != (b == length)) fail("cycle kind");
                Sequence qbar(s.rbegin() + (length - b), s.rbegin() + (length - a));
                for (int& v : qbar) v = -v;
                if (op.kind == 'f') {
                    std::reverse(s.begin() + b, s.end());
                    for (std::size_t c = b; c < length; ++c) s[c] = -s[c];
                }
                s.insert(s.end(), qbar.begin(), qbar.end());
            }
            cost.total += 1;
            if (!foldable(s)) fail("string not a foldable prefix");
        }
        if (s != target_) fail("history does not end in the target");
        if (!(cost == Cost{result.foldback_t, result.t_cycles, result.total})) {
            fail("cost");
        }
    }

    // The string a t-cycle appends: a listed source, or in block mode a
    // block T[x, x+y) inside one of the strings g or rc(g).
    Sequence appended(const Operation& op) const {
        if (!blocks_) return sources_.at(static_cast<std::size_t>(op.source));
        const std::uint32_t begin = static_cast<std::uint32_t>(op.x);
        const std::uint32_t end = begin + static_cast<std::uint32_t>(op.y);
        bool inside = false;
        for (const auto& [b, e] : strings_) inside |= b <= begin && begin < end && end <= e;
        if (!inside) fail("appended block is not a source");
        return Sequence(text_.begin() + begin, text_.begin() + end);
    }

    // Direct test of the definition of a foldable prefix: the tail after
    // the longest common prefix with t, or its reverse complement, occurs
    // in t.
    bool foldable(const Sequence& s) const {
        if (s.size() > n_) return false;
        auto phi = s.begin();
        auto at = target_.begin();
        while (phi != s.end() && *phi == *at) {
            ++phi;
            ++at;
        }
        if (phi == s.end()) return true;
        if (std::search(target_.begin(), target_.end(), phi, s.end()) != target_.end()) {
            return true;
        }
        auto equal_bar = [&](Sequence::const_iterator start) {
            auto back = s.end();
            for (auto it = start; back != phi; ++it) {
                --back;
                if (*it != -*back) return false;
            }
            return true;
        };
        const std::ptrdiff_t len = s.end() - phi;
        for (auto start = target_.begin(); target_.end() - start >= len; ++start) {
            if (equal_bar(start)) return true;
        }
        return false;
    }

    [[noreturn]] static void fail(const std::string& what) {
        throw std::logic_error("nbfb_simple replay failed: " + what);
    }

    // -------------------------------------------------------- state

    using HeapEntry = std::pair<Cost, std::uint32_t>;
    struct Later {
        bool operator()(const HeapEntry& a, const HeapEntry& b) const {
            return b.first < a.first;
        }
    };

    const std::vector<Sequence>& sources_;
    const Sequence& target_;
    std::uint32_t n_;
    const JoinTable* joins_;
    bool blocks_;
    Accounted accounted_;

    Sequence text_;
    // Block mode: the strings g_k and rc(g_k) as [begin, end) ranges of T,
    // and the index of every source line.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> strings_;
    std::unordered_map<Sequence, int, SequenceHash> listed_;
    std::size_t width_ = 0;
    std::vector<std::uint16_t> lce_;
    std::uint64_t base_ = 0;
    std::vector<std::uint64_t> prefix_;
    std::vector<std::uint64_t> power_;
    FlatMap windows_;
    struct TrieNode {
        std::vector<std::pair<int, std::uint32_t>> children;
        std::vector<std::uint32_t> sources;
    };
    std::vector<TrieNode> trie_;
    FlatMap root_;                       // first symbol -> trie child
    int symbols_ = 0;                    // largest |symbol| of t
    FlatMap bigrams_;                    // adjacent pairs of t and rc(t)
    FlatMap follow_index_;               // window -> index into follow_begin_
    std::vector<std::uint32_t> follow_begin_;
    std::vector<int> follow_symbols_;
    std::vector<std::uint32_t> stamp_;   // per symbol: follows the tail
    std::uint32_t epoch_ = 0;

    std::vector<Node> nodes_;
    FlatMap index_;
    std::vector<HeapEntry> heap_;  // binary heap ordered by Later
    std::vector<std::uint32_t> path_;
};

}  // namespace

Result nbfb_simple_dist(const PreparedSources& sources,
                        const std::vector<int>& target,
                        bool want_history,
                        MemoryBudget& budget,
                        const JoinTable* foldback_first,
                        bool blocks) {
    if (target.empty()) {
        Result result;
        result.reachable = true;
        return result;
    }
    SimpleSearch search(sources, target, budget, foldback_first, blocks);
    return search.solve(want_history);
}

}  // namespace nbfb
