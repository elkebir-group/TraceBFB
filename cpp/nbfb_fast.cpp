#include "nbfb_fast.hpp"

#include <algorithm>
#include <limits>
#include <utility>

namespace nbfb {

// ---------------------------------------------------------------- budget

MemoryBudget::MemoryBudget(std::size_t limit_bytes)
    : limit_(limit_bytes), used_(0), peak_(0) {}

void MemoryBudget::acquire(std::size_t bytes) {
    std::size_t now = used_.fetch_add(bytes, std::memory_order_relaxed) + bytes;
    if (now > limit_) {
        used_.fetch_sub(bytes, std::memory_order_relaxed);
        throw MemoryLimitExceeded(
            "memory limit exceeded: the search needs more than " +
            std::to_string(limit_) + " bytes");
    }
    std::size_t seen = peak_.load(std::memory_order_relaxed);
    while (now > seen &&
           !peak_.compare_exchange_weak(seen, now, std::memory_order_relaxed)) {
    }
}

void MemoryBudget::release(std::size_t bytes) {
    used_.fetch_sub(bytes, std::memory_order_relaxed);
}

PreparedSources::PreparedSources(std::vector<std::vector<int>> sources)
    : forward(std::move(sources)) {
    reverse.reserve(forward.size());
    for (const std::vector<int>& source : forward) {
        std::vector<int> other(source.rbegin(), source.rend());
        for (int& value : other) value = -value;
        reverse.push_back(std::move(other));
    }
}

namespace {

using i32 = std::int32_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

constexpr u32 NO_PARENT = std::numeric_limits<u32>::max();
constexpr u64 EMPTY_KEY = std::numeric_limits<u64>::max();

// The four packed fields of an occurrence or MINB key must stay below
// 2^16, and the source index must stay below 0xFFFF so that no key can
// collide with the empty-slot sentinel.
constexpr int MAX_INDEX = 65534;

inline u64 mix64(u64 value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

inline u64 pack4(int a, i32 b, i32 c, i32 d) {
    return (static_cast<u64>(static_cast<u32>(a)) << 48) |
           (static_cast<u64>(static_cast<u32>(b)) << 32) |
           (static_cast<u64>(static_cast<u32>(c)) << 16) |
           static_cast<u64>(static_cast<u32>(d));
}

// Open-addressing map from a packed 64-bit key to a pair of 32-bit
// values, used for the memoized OCC and MINB tables. Both tables are
// sparse (only the entries an expansion asks for are ever computed), so a
// hash map is the right structure; the dense index space would be
// O(m n^3).
class PairMap {
  public:
    explicit PairMap(MemoryBudget& budget) : budget_(&budget) {}
    PairMap(const PairMap&) = delete;
    PairMap& operator=(const PairMap&) = delete;
    ~PairMap() {
        if (bytes_ != 0) budget_->release(bytes_);
    }

    bool find(u64 key, i32& first, i32& second) const {
        if (slots_.empty()) return false;
        std::size_t position = static_cast<std::size_t>(mix64(key)) & mask_;
        while (true) {
            const Slot& slot = slots_[position];
            if (slot.key == EMPTY_KEY) return false;
            if (slot.key == key) {
                first = slot.first;
                second = slot.second;
                return true;
            }
            position = (position + 1) & mask_;
        }
    }

    void insert(u64 key, i32 first, i32 second) {
        if ((size_ + 1) * 2 > slots_.size()) grow();
        std::size_t position = static_cast<std::size_t>(mix64(key)) & mask_;
        while (slots_[position].key != EMPTY_KEY) {
            if (slots_[position].key == key) {
                slots_[position].first = first;
                slots_[position].second = second;
                return;
            }
            position = (position + 1) & mask_;
        }
        slots_[position] = Slot{key, first, second};
        ++size_;
    }

  private:
    struct Slot {
        u64 key;
        i32 first;
        i32 second;
    };

    void grow() {
        std::size_t capacity = slots_.empty() ? 1024 : slots_.size() * 2;
        std::size_t wanted = capacity * sizeof(Slot);
        budget_->acquire(wanted);
        std::vector<Slot> next(capacity, Slot{EMPTY_KEY, 0, 0});
        std::size_t next_mask = capacity - 1;
        for (const Slot& slot : slots_) {
            if (slot.key == EMPTY_KEY) continue;
            std::size_t position =
                static_cast<std::size_t>(mix64(slot.key)) & next_mask;
            while (next[position].key != EMPTY_KEY) {
                position = (position + 1) & next_mask;
            }
            next[position] = slot;
        }
        if (bytes_ != 0) budget_->release(bytes_);
        bytes_ = wanted;
        slots_.swap(next);
        mask_ = next_mask;
    }

    MemoryBudget* budget_;
    std::vector<Slot> slots_;
    std::size_t mask_ = 0;
    std::size_t size_ = 0;
    std::size_t bytes_ = 0;
};

// Interval union-find over lo..hi supporting "smallest unvisited value at
// or above v" and "mark v visited" in amortized near-constant time, the
// C++ form of the Python _NextUnvisited.
inline i32 nu_next(i32* parent, i32 lo, i32 hi, i32 value) {
    if (value < lo) value = lo;
    if (value > hi) return hi + 1;
    i32 root = value;
    while (parent[root - lo] != root) root = parent[root - lo];
    while (parent[value - lo] != root) {
        i32 next = parent[value - lo];
        parent[value - lo] = root;
        value = next;
    }
    return root;
}

inline void nu_remove(i32* parent, i32 lo, i32 value) {
    parent[value - lo] = value + 1;
}

inline bool nu_unvisited(const i32* parent, i32 lo, i32 value) {
    return parent[value - lo] == value;
}

constexpr std::uint8_t KIND_E = 0;
constexpr std::uint8_t KIND_W = 1;
constexpr std::uint8_t KIND_T = 2;

// One visited index tuple. The vector of records doubles as the
// breadth-first queue (records are never removed) and as the parent-
// pointer table used to reconstruct one optimal history.
struct Record {
    u32 parent;
    u32 dist;
    i32 phi;  // E and W: phi; T: phi0
    i32 a;    // W: window start; T: source index
    i32 b;    // W: window end
    i32 op_x;
    i32 op_y;
    std::uint16_t mask;  // state: foldback_t * (n + 1) + t-BFB cycles so far
    std::uint8_t kind;
    char op_kind;
};

// The per-layer (per t-BFB count) visited structures. Every table is
// allocated the first time the layer touches a node of that kind, and
// the window blocks are allocated one (phi, b) pair at a time, so an
// instance only pays for the tuples the search actually reaches.
struct Layer {
    std::vector<i32> e_parent;   // values 0..n+1
    std::vector<i32> w_index;    // (n+1)*(n+1) block ids, -1 when absent
    std::vector<i32> w_pool;     // blocks of n+1 entries, values 1..n+1
    std::vector<std::uint8_t> t_seen;  // (n+1)*m flags
    std::size_t bytes = 0;
};

class Search {
  public:
    Search(const PreparedSources& sources, const std::vector<int>& target,
           MemoryBudget& budget, const JoinTable* joins = nullptr);
    ~Search();

    Result solve(bool want_history);

  private:
    // -- allocation ----------------------------------------------------
    template <class T>
    void alloc(std::vector<T>& vector, std::size_t count, T value) {
        std::size_t bytes = count * sizeof(T);
        budget_->acquire(bytes);
        owned_ += bytes;
        vector.assign(count, value);
    }

    // -- tables --------------------------------------------------------
    void build_tables();
    std::vector<i32> lce_table(const std::vector<int>& y,
                               const std::vector<int>& z);
    void build_lce_tx(int index);
    void build_keep_t(int index);

    void occ(int zkey, i32 i, i32 lambda, i32 s, i32& best, i32& pos);
    void minb(int zkey, i32 phi, i32 s, i32 mu, i32& answer, i32& pos);

    // -- layers --------------------------------------------------------
    Layer& layer(int mask);
    i32* e_parent(Layer& target_layer);
    i32* w_parent(Layer& target_layer, i32 phi, i32 b);
    std::uint8_t* t_seen(Layer& target_layer);

    // -- search --------------------------------------------------------
    void reset(int max_translocations, int max_foldback = 0, bool count_t = true);
    long long run();
    std::vector<Operation> history_of(long long goal) const;

    // The search state packs the number f of t-BFB cycles that created a
    // fold-back-shaped join and the number t of t-BFB cycles as
    // f * (n + 1) + t; without a JoinTable f stays 0 and the state is t.
    int t_of(int state) const { return state % state_stride_; }
    int f_of(int state) const { return state / state_stride_; }
    bool may_append(int state) const {
        return t_of(state) < max_translocations_;
    }
    // State after a t-BFB cycle appending a source that starts with v to a
    // sequence ending with u (u = 0 for the empty sequence), or -1 when the
    // bound on fold-back-shaped joins forbids it.
    int after_t(int state, int u, int v) const {
        const int fb = (joins_ != nullptr && u != 0 && joins_->is_foldback(u, v)) ? 1 : 0;
        if (f_of(state) + fb > max_foldback_) return -1;
        return state + (count_t_ ? 1 : 0) + fb * state_stride_;
    }
    void push(std::uint8_t kind, int mask, i32 phi, i32 a, i32 b, u32 parent,
              char op_kind, i32 op_x, i32 op_y);
    void visit_E(int mask, i32 phi, u32 parent, char op_kind, i32 op_x,
                 i32 op_y);
    void visit_W(int mask, i32 phi, i32 a, i32 b, u32 parent, char op_kind,
                 i32 op_x, i32 op_y);
    void visit_T(int mask, i32 phi0, i32 index, u32 parent, char op_kind,
                 i32 op_x, i32 op_y);
    void expand_empty(int mask, i32 phi, u32 key);
    void expand_tail(const Record& record, u32 key);

    const PreparedSources* sources_;
    MemoryBudget* budget_;
    std::size_t owned_ = 0;

    i32 n_ = 0;
    int m_ = 0;
    std::vector<int> t_;     // target, 0-based storage, 1-based indices
    std::vector<int> tbar_;  // reverse complement of the target
    std::vector<i32> src_len_;

    std::vector<i32> lce_tt_;  // (n+2) x (n+2)
    std::vector<i32> lce_bt_;  // (n+2) x (n+2)
    std::vector<i32> lce_tb_;  // (n+2) x (n+2)
    std::size_t stride_n_ = 0;

    // Row 1 of LCE(x, t) and column 1 of LCE(t, xbar): the only slices of
    // those tables the algorithm reads.
    std::vector<i32> lce_xt1_;     // m x (n+2)
    std::vector<i32> lce_txbar1_;  // m x (n+2)

    // Full LCE(t, x) tables, built on demand for the sources that occur
    // as source nodes.
    std::vector<std::vector<i32>> lce_tx_;

    std::vector<i32> xbar_off_;  // CSR of the occurrences of every xbar
    std::vector<i32> xbar_pos_;
    std::vector<int> foldable_;

    struct Entry {
        i32 length;
        i32 pos;
        i32 source;
    };
    std::vector<i32> keep_w_off_;  // CSR over window starts a = 1..n
    std::vector<Entry> keep_w_;
    std::vector<std::vector<Entry>> keep_t_;  // built on demand
    std::vector<char> keep_t_ready_;

    PairMap occ_map_;
    PairMap minb_map_;

    std::vector<Record> records_;
    std::size_t record_bytes_ = 0;
    std::vector<Layer> layers_;
    int max_translocations_ = 0;
    const JoinTable* joins_ = nullptr;
    int max_foldback_ = 0;
    // false only in the feasibility runs of the fold-back-first objective: t-BFB
    // cycles are then not counted, so a state is (f, tuple) and reachability with at
    // most f fold-back-shaped t-joins is decided in one layer per f
    bool count_t_ = true;
    int state_stride_ = 1;
};

Search::Search(const PreparedSources& sources, const std::vector<int>& target,
               MemoryBudget& budget, const JoinTable* joins)
    : sources_(&sources),
      budget_(&budget),
      occ_map_(budget),
      minb_map_(budget),
      joins_(joins) {
    n_ = static_cast<i32>(target.size());
    state_stride_ = static_cast<int>(n_) + 1;
    m_ = static_cast<int>(sources.forward.size());
    if (n_ > MAX_INDEX || m_ > MAX_INDEX) {
        throw std::invalid_argument(
            "instance too large: the target and the number of sources must "
            "both stay below 65535");
    }
    t_ = target;
    tbar_.resize(target.size());
    for (std::size_t i = 0; i < target.size(); ++i) {
        tbar_[i] = -target[target.size() - 1 - i];
    }
    try {
        build_tables();
    } catch (...) {
        // The destructor does not run when a constructor throws, so the
        // budget is returned here before the exception leaves.
        budget_->release(owned_);
        owned_ = 0;
        throw;
    }
}

Search::~Search() {
    for (Layer& item : layers_) budget_->release(item.bytes);
    budget_->release(record_bytes_);
    budget_->release(owned_);
}

// ---------------------------------------------------------------- tables

std::vector<i32> Search::lce_table(const std::vector<int>& y,
                                   const std::vector<int>& z) {
    const std::size_t ly = y.size();
    const std::size_t lz = z.size();
    const std::size_t stride = lz + 2;
    std::vector<i32> table;
    alloc(table, (ly + 2) * stride, i32(0));
    for (std::size_t i = ly; i >= 1; --i) {
        const int yi = y[i - 1];
        i32* row = table.data() + i * stride;
        const i32* next = table.data() + (i + 1) * stride;
        for (std::size_t j = lz; j >= 1; --j) {
            if (yi == z[j - 1]) row[j] = 1 + next[j + 1];
        }
    }
    return table;
}

void Search::build_tables() {
    stride_n_ = static_cast<std::size_t>(n_) + 2;
    lce_tt_ = lce_table(t_, t_);
    lce_bt_ = lce_table(tbar_, t_);
    lce_tb_ = lce_table(t_, tbar_);

    src_len_.resize(m_);
    for (int index = 0; index < m_; ++index) {
        src_len_[index] = static_cast<i32>(sources_->forward[index].size());
        if (src_len_[index] == 0) {
            throw std::invalid_argument("source " + std::to_string(index) +
                                        " is empty");
        }
    }

    alloc(lce_xt1_, static_cast<std::size_t>(m_) * stride_n_, i32(0));
    alloc(lce_txbar1_, static_cast<std::size_t>(m_) * stride_n_, i32(0));
    for (int index = 0; index < m_; ++index) {
        const std::vector<int>& x = sources_->forward[index];
        const std::vector<int>& xbar = sources_->reverse[index];
        const i32 length = src_len_[index];
        i32* xt = lce_xt1_.data() + static_cast<std::size_t>(index) * stride_n_;
        i32* txbar =
            lce_txbar1_.data() + static_cast<std::size_t>(index) * stride_n_;
        for (i32 j = 1; j <= n_; ++j) {
            i32 k = 0;
            while (k < length && j + k <= n_ && x[k] == t_[j + k - 1]) ++k;
            xt[j] = k;
            k = 0;
            while (k < length && j + k <= n_ && t_[j + k - 1] == xbar[k]) ++k;
            txbar[j] = k;
        }
    }

    // Occurrences of every reverse-complemented source in the target.
    alloc(xbar_off_, static_cast<std::size_t>(m_) + 1, i32(0));
    std::vector<i32> positions;
    for (int index = 0; index < m_; ++index) {
        const i32 length = src_len_[index];
        const i32* txbar =
            lce_txbar1_.data() + static_cast<std::size_t>(index) * stride_n_;
        for (i32 p = 1; p <= n_ - length + 1; ++p) {
            if (txbar[p] >= length) positions.push_back(p);
        }
        xbar_off_[index + 1] = static_cast<i32>(positions.size());
    }
    alloc(xbar_pos_, positions.size(), i32(0));
    std::copy(positions.begin(), positions.end(), xbar_pos_.begin());

    for (int y = 0; y < m_; ++y) {
        if (xbar_off_[y + 1] > xbar_off_[y]) foldable_.push_back(y);
    }
    budget_->acquire(foldable_.capacity() * sizeof(int));
    owned_ += foldable_.capacity() * sizeof(int);

    // For every window start a, the sources that can be appended while
    // keeping the matched prefix, by decreasing extension length.
    const std::size_t foldable_count = foldable_.size();
    std::vector<i32> best_length;
    std::vector<i32> best_pos;
    alloc(best_length, foldable_count * stride_n_, i32(0));
    alloc(best_pos, foldable_count * stride_n_, i32(0));
    for (std::size_t rank = 0; rank < foldable_count; ++rank) {
        const int y = foldable_[rank];
        const i32 ly = src_len_[y];
        i32* length_row = best_length.data() + rank * stride_n_;
        i32* pos_row = best_pos.data() + rank * stride_n_;
        for (i32 index = xbar_off_[y]; index < xbar_off_[y + 1]; ++index) {
            const i32 p = xbar_pos_[index];
            const i32* row =
                lce_tt_.data() + static_cast<std::size_t>(p + ly) * stride_n_;
            for (i32 a = 1; a <= n_; ++a) {
                if (row[a] > length_row[a]) {
                    length_row[a] = row[a];
                    pos_row[a] = p;
                }
            }
        }
    }
    alloc(keep_w_off_, static_cast<std::size_t>(n_) + 2, i32(0));
    std::vector<Entry> entries;
    std::vector<Entry> row;
    for (i32 a = 1; a <= n_; ++a) {
        row.clear();
        for (std::size_t rank = 0; rank < foldable_count; ++rank) {
            const i32 length = best_length[rank * stride_n_ + a];
            if (length >= 1) {
                row.push_back(
                    Entry{length, best_pos[rank * stride_n_ + a],
                          static_cast<i32>(foldable_[rank])});
            }
        }
        std::stable_sort(row.begin(), row.end(),
                         [](const Entry& left, const Entry& right) {
                             return left.length > right.length;
                         });
        entries.insert(entries.end(), row.begin(), row.end());
        keep_w_off_[a + 1] = static_cast<i32>(entries.size());
    }
    alloc(keep_w_, entries.size(), Entry{0, 0, 0});
    std::copy(entries.begin(), entries.end(), keep_w_.begin());
    budget_->release(best_length.size() * sizeof(i32));
    budget_->release(best_pos.size() * sizeof(i32));
    owned_ -= best_length.size() * sizeof(i32);
    owned_ -= best_pos.size() * sizeof(i32);
    best_length.clear();
    best_length.shrink_to_fit();
    best_pos.clear();
    best_pos.shrink_to_fit();

    lce_tx_.resize(m_);
    keep_t_.resize(m_);
    keep_t_ready_.assign(m_, 0);
}

void Search::build_lce_tx(int index) {
    if (!lce_tx_[index].empty()) return;
    lce_tx_[index] = lce_table(t_, sources_->forward[index]);
}

void Search::build_keep_t(int index) {
    if (keep_t_ready_[index]) return;
    keep_t_ready_[index] = 1;
    const i32* txbar =
        lce_txbar1_.data() + static_cast<std::size_t>(index) * stride_n_;
    std::vector<Entry>& row = keep_t_[index];
    for (int y : foldable_) {
        const i32 ly = src_len_[y];
        i32 best = 0;
        i32 best_pos = 0;
        for (i32 slot = xbar_off_[y]; slot < xbar_off_[y + 1]; ++slot) {
            const i32 p = xbar_pos_[slot];
            const i32 value = txbar[p + ly];
            if (value > best) {
                best = value;
                best_pos = p;
            }
        }
        if (best >= 1) row.push_back(Entry{best, best_pos, static_cast<i32>(y)});
    }
    std::stable_sort(row.begin(), row.end(),
                     [](const Entry& left, const Entry& right) {
                         return left.length > right.length;
                     });
    std::size_t bytes = row.capacity() * sizeof(Entry);
    budget_->acquire(bytes);
    owned_ += bytes;
}

void Search::occ(int zkey, i32 i, i32 lambda, i32 s, i32& best, i32& pos) {
    const u64 key = pack4(zkey, i, lambda, s);
    if (occ_map_.find(key, best, pos)) return;

    const i32* zt;
    std::size_t stride;
    if (zkey == m_) {
        zt = lce_tb_.data();
        stride = stride_n_;
    } else {
        build_lce_tx(zkey);
        zt = lce_tx_[zkey].data();
        stride = static_cast<std::size_t>(src_len_[zkey]) + 2;
    }
    best = -1;
    pos = i;
    const i32* tt = lce_tt_.data();
    const i32 last = n_ - lambda + 1;
    for (i32 p = 1; p <= last; ++p) {
        if (tt[static_cast<std::size_t>(p) * stride_n_ + i] >= lambda) {
            const i32 value =
                zt[static_cast<std::size_t>(p + lambda) * stride + s];
            if (value > best) {
                best = value;
                pos = p;
            }
        }
    }
    occ_map_.insert(key, best, pos);
}

void Search::minb(int zkey, i32 phi, i32 s, i32 mu, i32& answer, i32& pos) {
    const u64 key = pack4(zkey, phi, s, mu);
    if (minb_map_.find(key, answer, pos)) return;

    // OCC_z(beta + 1, phi - beta, s) is nondecreasing in beta.
    i32 lo = 0;
    i32 hi = phi - 1;
    answer = phi;
    pos = -1;
    while (lo <= hi) {
        const i32 mid = (lo + hi) / 2;
        i32 value;
        i32 position;
        occ(zkey, mid + 1, phi - mid, s, value, position);
        if (value >= mu) {
            answer = mid;
            pos = position;
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }
    minb_map_.insert(key, answer, pos);
}

// ---------------------------------------------------------------- layers

Layer& Search::layer(int mask) {
    if (static_cast<std::size_t>(mask) >= layers_.size()) {
        layers_.resize(static_cast<std::size_t>(mask) + 1);
    }
    return layers_[mask];
}

i32* Search::e_parent(Layer& target_layer) {
    if (target_layer.e_parent.empty()) {
        const std::size_t count = static_cast<std::size_t>(n_) + 2;
        budget_->acquire(count * sizeof(i32));
        target_layer.bytes += count * sizeof(i32);
        target_layer.e_parent.resize(count);
        for (std::size_t v = 0; v < count; ++v) {
            target_layer.e_parent[v] = static_cast<i32>(v);
        }
    }
    return target_layer.e_parent.data();
}

i32* Search::w_parent(Layer& target_layer, i32 phi, i32 b) {
    // The successor rules never leave these ranges; the check turns a
    // hypothetical rule error into an exception rather than a silent
    // out-of-bounds write, at the cost of two comparisons per window node.
    if (phi < 0 || phi > n_ || b < 1 || b > n_) {
        throw std::logic_error("window node W(" + std::to_string(phi) + ", *, " +
                               std::to_string(b) + ") is out of range");
    }
    const std::size_t width = static_cast<std::size_t>(n_) + 1;
    if (target_layer.w_index.empty()) {
        const std::size_t count = width * width;
        budget_->acquire(count * sizeof(i32));
        target_layer.bytes += count * sizeof(i32);
        target_layer.w_index.assign(count, -1);
    }
    const std::size_t slot =
        static_cast<std::size_t>(phi) * width + static_cast<std::size_t>(b - 1);
    i32 block = target_layer.w_index[slot];
    if (block < 0) {
        const std::size_t block_size = width;  // values 1..n+1
        budget_->acquire(block_size * sizeof(i32));
        target_layer.bytes += block_size * sizeof(i32);
        block = static_cast<i32>(target_layer.w_pool.size() / block_size);
        target_layer.w_pool.resize(target_layer.w_pool.size() + block_size);
        i32* parent = target_layer.w_pool.data() +
                      static_cast<std::size_t>(block) * block_size;
        for (std::size_t v = 0; v < block_size; ++v) {
            parent[v] = static_cast<i32>(v) + 1;
        }
        target_layer.w_index[slot] = block;
    }
    return target_layer.w_pool.data() +
           static_cast<std::size_t>(block) * width;
}

std::uint8_t* Search::t_seen(Layer& target_layer) {
    if (target_layer.t_seen.empty()) {
        const std::size_t count =
            (static_cast<std::size_t>(n_) + 1) * static_cast<std::size_t>(m_);
        budget_->acquire(count);
        target_layer.bytes += count;
        target_layer.t_seen.assign(count, 0);
    }
    return target_layer.t_seen.data();
}

// ---------------------------------------------------------------- search

void Search::reset(int max_translocations, int max_foldback, bool count_t) {
    max_translocations_ = max_translocations;
    max_foldback_ = max_foldback;
    count_t_ = count_t;
    for (Layer& item : layers_) budget_->release(item.bytes);
    layers_.clear();
    layers_.shrink_to_fit();
    budget_->release(record_bytes_);
    record_bytes_ = 0;
    records_.clear();
    records_.shrink_to_fit();
}

void Search::push(std::uint8_t kind, int mask, i32 phi, i32 a, i32 b,
                  u32 parent, char op_kind, i32 op_x, i32 op_y) {
    if (records_.size() == records_.capacity()) {
        const std::size_t next =
            records_.capacity() == 0 ? 1024 : records_.capacity() * 2;
        const std::size_t bytes = next * sizeof(Record);
        budget_->acquire(bytes);
        budget_->release(record_bytes_);
        record_bytes_ = bytes;
        records_.reserve(next);
    }
    Record record;
    record.parent = parent;
    record.dist = parent == NO_PARENT ? 0u : records_[parent].dist + 1u;
    record.phi = phi;
    record.a = a;
    record.b = b;
    record.op_x = op_x;
    record.op_y = op_y;
    record.mask = static_cast<std::uint16_t>(mask);
    record.kind = kind;
    record.op_kind = op_kind;
    records_.push_back(record);
}

void Search::visit_E(int mask, i32 phi, u32 parent, char op_kind, i32 op_x,
                     i32 op_y) {
    if (phi < 0 || phi > n_) {
        throw std::logic_error("empty-tail node E(" + std::to_string(phi) +
                               ") is out of range");
    }
    i32* family = e_parent(layer(mask));
    if (!nu_unvisited(family, 0, phi)) return;
    nu_remove(family, 0, phi);
    push(KIND_E, mask, phi, 0, 0, parent, op_kind, op_x, op_y);
}

void Search::visit_W(int mask, i32 phi, i32 a, i32 b, u32 parent, char op_kind,
                     i32 op_x, i32 op_y) {
    i32* family = w_parent(layer(mask), phi, b);
    if (!nu_unvisited(family, 1, a)) return;
    nu_remove(family, 1, a);
    push(KIND_W, mask, phi, a, b, parent, op_kind, op_x, op_y);
}

void Search::visit_T(int mask, i32 phi0, i32 index, u32 parent, char op_kind,
                     i32 op_x, i32 op_y) {
    if (phi0 < 0 || phi0 > n_ || index < 0 || index >= m_) {
        throw std::logic_error("source node T(" + std::to_string(phi0) + ", " +
                               std::to_string(index) + ") is out of range");
    }
    std::uint8_t* seen = t_seen(layer(mask));
    const std::size_t slot = static_cast<std::size_t>(phi0) *
                                 static_cast<std::size_t>(m_) +
                             static_cast<std::size_t>(index);
    if (seen[slot] != 0) return;
    seen[slot] = 1;
    push(KIND_T, mask, phi0, index, 0, parent, op_kind, op_x, op_y);
}

void Search::expand_empty(int mask, i32 phi, u32 key) {
    if (may_append(mask)) {
        const int last = phi >= 1 ? t_[phi - 1] : 0;
        for (int index = 0; index < m_; ++index) {
            const int next = after_t(mask, last, sources_->forward[index][0]);
            if (next < 0) continue;
            const i32 length = src_len_[index];
            const i32 lce =
                lce_xt1_[static_cast<std::size_t>(index) * stride_n_ + phi + 1];
            const i32 j = std::min(length, lce);
            if (j >= 1) {
                if (j == length) {
                    visit_E(next, phi + length, key, 't', index, -1);
                } else {
                    visit_T(next, phi, index, key, 't', index, -1);
                }
            } else if (xbar_off_[index + 1] > xbar_off_[index]) {
                // The prefix is kept; the tail x is foldable as xbar occurs.
                const i32 p = xbar_pos_[xbar_off_[index]];
                visit_W(next, phi, p, p + length - 1, key, 't', index, -1);
            }
        }
    }

    const i32* bt = lce_bt_.data();
    const i32* bt_prefix =
        bt + static_cast<std::size_t>(n_ - phi + 1) * stride_n_;
    for (i32 beta = 0; beta < phi; ++beta) {
        // p-BFB: q = t[beta+1..phi].
        {
            const i32 M = std::min(phi - beta, bt_prefix[phi + 1]);
            if (M >= 1) {
                if (M == phi - beta) {
                    visit_E(mask, phi + M, key, 'p', beta, phi);
                } else {
                    visit_W(mask, phi + M, beta + 1, phi - M, key, 'p', beta,
                            phi);
                }
            }
        }
        for (i32 gamma = beta + 1; gamma < phi; ++gamma) {
            // i-BFB: q = t[beta+1..gamma], r = t[gamma+1..phi].
            {
                const i32 M = std::min(
                    gamma - beta,
                    bt[static_cast<std::size_t>(n_ - gamma + 1) * stride_n_ +
                       phi + 1]);
                if (M >= 1) {
                    if (M == gamma - beta) {
                        visit_E(mask, phi + M, key, 'i', beta, gamma);
                    } else {
                        visit_W(mask, phi + M, beta + 1, gamma - M, key, 'i',
                                beta, gamma);
                    }
                }
            }
            // f-BFB: q = t[beta+1..c], r = t[c+1..phi], c = gamma.
            {
                const i32 c = gamma;
                const i32 M = std::min(phi - beta, bt_prefix[c + 1]);
                if (c + M > phi) {
                    if (M == phi - beta) {
                        visit_E(mask, c + M, key, 'f', beta, c);
                    } else {
                        visit_W(mask, c + M, beta + 1, phi - M, key, 'f', beta,
                                c);
                    }
                }
            }
        }
    }
}

void Search::expand_tail(const Record& record, u32 key) {
    const int mask = record.mask;
    i32 phi;
    i32 lambda;
    // The number of symbols of the reverse-complemented tail matched by
    // t[c+1..n] is base_c[(c + 1) * step_c]: a column of LCE(t, t) for a
    // window node and a row of LCE(t, xbar) for a source node.
    const i32* base_c;
    std::size_t step_c;
    int zkey;
    i32 s;
    i32 a = 0;
    int index = 0;

    if (record.kind == KIND_W) {
        phi = record.phi;
        a = record.a;
        lambda = record.b - record.a + 1;
        base_c = lce_tt_.data() + a;
        step_c = stride_n_;
        zkey = m_;
        s = n_ - record.b + 1;
    } else {
        index = record.a;
        const i32 length = src_len_[index];
        const i32 j = std::min(
            length,
            lce_xt1_[static_cast<std::size_t>(index) * stride_n_ + record.phi +
                     1]);
        phi = record.phi + j;
        lambda = length - j;
        base_c = lce_txbar1_.data() + static_cast<std::size_t>(index) * stride_n_;
        step_c = 1;
        zkey = index;
        s = j + 1;
    }

    const i32* bt_prefix =
        lce_bt_.data() + static_cast<std::size_t>(n_ - phi + 1) * stride_n_;

    for (i32 c = 1; c <= phi; ++c) {
        const i32 extension = base_c[static_cast<std::size_t>(c + 1) * step_c];
        const i32 M = std::min(lambda, extension);
        if (M == 0) continue;

        if (M == lambda) {
            const i32 L = bt_prefix[c + lambda + 1];
            // Full absorption, empty-tail successors E(phi') with
            // beta = phi - (phi' - c - lambda).
            const i32 lo = lambda + phi + 1;
            const i32 hi = c + lambda + std::min(phi, L);
            if (lo <= hi) {
                i32* family = e_parent(layer(mask));
                i32 v = nu_next(family, 0, n_, lo);
                while (v <= hi) {
                    nu_remove(family, 0, v);
                    push(KIND_E, mask, v, 0, 0, key, 'f',
                         phi - (v - c - lambda), c);
                    v = nu_next(family, 0, n_, v + 1);
                }
            }
            // Full absorption, window successors W(c+lambda+L, a', phi-L)
            // with beta = a' - 1.
            if (c + lambda + L > phi) {
                const i32 hi_a = std::min(c, phi - L);
                if (hi_a >= 1) {
                    const i32 new_phi = c + lambda + L;
                    const i32 new_b = phi - L;
                    // No further window block is allocated inside the loop,
                    // so the block pointer stays valid.
                    i32* family = w_parent(layer(mask), new_phi, new_b);
                    i32 v = nu_next(family, 1, n_, 1);
                    while (v <= hi_a) {
                        nu_remove(family, 1, v);
                        push(KIND_W, mask, new_phi, v, new_b, key, 'f', v - 1,
                             c);
                        v = nu_next(family, 1, n_, v + 1);
                    }
                }
            }
            continue;
        }

        // Partial absorption.
        if (c + M <= phi) continue;
        const i32 mu = lambda - M;
        i32 beta_min;
        i32 P;
        minb(zkey, phi, s, mu, beta_min, P);
        if (beta_min > c - 1) continue;
        const i32 new_phi = c + M;
        const i32 new_b = P - beta_min + phi + mu - 1;
        const i32 hi_a = P + c - 1 - beta_min;
        i32* family = w_parent(layer(mask), new_phi, new_b);
        i32 v = nu_next(family, 1, n_, P);
        while (v <= hi_a) {
            nu_remove(family, 1, v);
            push(KIND_W, mask, new_phi, v, new_b, key, 'f', beta_min + (v - P),
                 c);
            v = nu_next(family, 1, n_, v + 1);
        }
    }

    // t-BFB cycles that keep the matched prefix: appending y to the tail u
    // is admitted when ybar ubar occurs in t, and the successor is the
    // window node whose window spells ybar ubar.
    if (may_append(mask)) {
        // Last segment of the tail: a window spells the reverse complement
        // of the tail, and a source node's tail is a suffix of its source.
        const int tail_last = record.kind == KIND_W
                                  ? -t_[a - 1]
                                  : sources_->forward[index].back();
        const Entry* first;
        const Entry* last;
        if (record.kind == KIND_W) {
            first = keep_w_.data() + keep_w_off_[a];
            last = keep_w_.data() + keep_w_off_[a + 1];
        } else {
            build_keep_t(index);
            const std::vector<Entry>& row = keep_t_[index];
            first = row.data();
            last = row.data() + row.size();
        }
        for (const Entry* entry = first; entry != last; ++entry) {
            if (entry->length < lambda) break;
            const int next =
                after_t(mask, tail_last, sources_->forward[entry->source][0]);
            if (next < 0) continue;
            visit_W(next, phi, entry->pos,
                    entry->pos + src_len_[entry->source] + lambda - 1, key, 't',
                    entry->source, -1);
        }
    }
}

long long Search::run() {
    // The start node E(0) with no t-BFB cycle applied.
    i32* family = e_parent(layer(0));
    nu_remove(family, 0, 0);
    push(KIND_E, 0, 0, 0, 0, NO_PARENT, 'p', -1, -1);

    for (std::size_t head = 0; head < records_.size(); ++head) {
        const std::uint8_t kind = records_[head].kind;
        if (kind == KIND_E) {
            if (records_[head].phi == n_) return static_cast<long long>(head);
            expand_empty(records_[head].mask, records_[head].phi,
                         static_cast<u32>(head));
        } else {
            // expand_tail may reallocate records_, so it works on a copy.
            const Record snapshot = records_[head];
            expand_tail(snapshot, static_cast<u32>(head));
        }
    }
    return -1;
}

std::vector<Operation> Search::history_of(long long goal) const {
    std::vector<Operation> history;
    u32 node = static_cast<u32>(goal);
    while (records_[node].parent != NO_PARENT) {
        const Record& current = records_[node];
        Operation operation;
        operation.kind = current.op_kind;
        if (current.op_kind == 't') {
            operation.source = current.op_x;
            operation.foldback =
                f_of(current.mask) > f_of(records_[current.parent].mask);
        } else {
            operation.x = current.op_x;
            operation.y = current.op_y;
        }
        history.push_back(operation);
        node = current.parent;
    }
    std::reverse(history.begin(), history.end());
    return history;
}

Result Search::solve(bool want_history) {
    long long goal = -1;
    Result result;
    if (joins_ == nullptr) {
        for (int bound = 1; bound <= n_; ++bound) {
            reset(bound);
            goal = run();
            if (goal >= 0) break;
        }
    } else {
        // Deepen the bound on fold-back-shaped t-BFB joins; for the first
        // feasible bound, deepen the bound on t-BFB cycles.
        const int cap = std::min<int>(n_, 65535 / state_stride_ - 1);
        for (int fb = 0; fb <= cap && goal < 0; ++fb) {
            reset(n_, fb, false);  // is the target reachable with at most fb such joins?
            const long long any = run();
            if (any < 0) continue;
            // With the t-BFB count unbounded, run() is a breadth-first
            // search, so this is the shortest history with at most fb
            // fold-back-shaped joins.
            const int shortest = static_cast<int>(records_[any].dist);
            result.fb_total = shortest;
            result.fb_total_t = -1;
            for (int bound = 1; bound <= n_; ++bound) {
                reset(bound, fb);
                const long long found = run();
                if (found < 0) continue;
                if (goal < 0) {
                    goal = found;
                    result.t_cycles = t_of(records_[goal].mask);
                    result.foldback_t = f_of(records_[goal].mask);
                    result.total = static_cast<int>(records_[goal].dist);
                    if (want_history) result.history = history_of(goal);
                }
                // The fewest t-BFB cycles among the shortest histories.
                if (static_cast<int>(records_[found].dist) == shortest) {
                    result.fb_total_t = t_of(records_[found].mask);
                    if (want_history) result.fb_total_history = history_of(found);
                    break;
                }
            }
        }
        if (goal < 0 && cap < n_) {
            throw std::invalid_argument(
                "target too long for the fold-back-first objective");
        }
        if (goal >= 0) {
            result.reachable = true;
            return result;
        }
    }
    if (goal < 0) {
        result.reachable = false;
        return result;
    }
    const Record& end = records_[static_cast<std::size_t>(goal)];
    result.reachable = true;
    result.t_cycles = t_of(end.mask);
    result.foldback_t = f_of(end.mask);
    result.total = static_cast<int>(end.dist);
    if (want_history) result.history = history_of(goal);
    return result;
}

}  // namespace

Result nbfb_fast_dist(const PreparedSources& sources,
                      const std::vector<int>& target, bool want_history,
                      MemoryBudget& budget, const JoinTable* foldback_first) {
    for (std::size_t index = 0; index < sources.forward.size(); ++index) {
        if (sources.forward[index].empty()) {
            throw std::invalid_argument("source " + std::to_string(index) +
                                        " is empty");
        }
    }
    if (target.empty()) {
        Result result;
        result.reachable = true;
        result.t_cycles = 0;
        result.total = 0;
        return result;
    }
    Search search(sources, target, budget, foldback_first);
    return search.solve(want_history);
}

Result nbfb_fast_dist(const std::vector<std::vector<int>>& sources,
                      const std::vector<int>& target, bool want_history) {
    PreparedSources prepared(sources);
    MemoryBudget budget(std::numeric_limits<std::size_t>::max());
    return nbfb_fast_dist(prepared, target, want_history, budget);
}

}  // namespace nbfb
