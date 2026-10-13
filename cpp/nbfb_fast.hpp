// Nested BFB distance over a reusable source set: C++17 port of
// src/bfb_dist/nbfb_fast.py (NBFB_fast_dist).
//
// A sequence of oriented genomic segments is represented by a vector of
// signed segment identifiers: every distinct segment interval gets a
// positive identifier, and the reverse complement of a segment is the
// negated identifier (a whole sequence is reverse complemented by
// reversing the order and negating every identifier).
//
// The search is the breadth-first search of the paper appendix "An
// O(2^m n^4)-Time Algorithm for NBFB Distance", over the index tuples
// E(phi), W(phi, a, b) and T(phi0, x), with the state extended by the
// number of t-BFB cycles applied so far. The driver deepens a bound on
// that number, so the returned distance is minimum first in t-BFB cycles
// and then in the total number of cycles.

#ifndef NBFB_FAST_HPP
#define NBFB_FAST_HPP

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace nbfb {

// One BFB cycle of a history, with the parameters of the Python
// bfb_dist.history.Operation. For kinds 'p', 'i' and 'f' the sequence the
// cycle acts on decomposes as P = v[:x], Q = v[x:y], R = v[y:]; for kind
// 't' the field source is the index of the appended source.
struct Operation {
    char kind = 'p';
    int x = -1;
    int y = -1;
    int source = -1;
    bool foldback = false;  // t-BFB cycle whose join is fold-back-shaped
};

struct Result {
    bool reachable = false;  // false means the distance is infinite
    int foldback_t = 0;      // t-BFB cycles creating a fold-back-shaped join
    int t_cycles = 0;        // number of t-BFB cycles of the optimum
    int total = 0;           // total number of cycles of the optimum
    std::vector<Operation> history;  // filled only when requested

    // Objective (fold-back-shaped t-BFB cycles, total cycles, t-BFB cycles):
    // the shortest history with the minimum number of fold-back-shaped
    // joins, and the fewest t-BFB cycles among those. Filled only when a
    // JoinTable is passed; equals (foldback_t, total, t_cycles) unless a
    // history with more t-BFB cycles is shorter.
    int fb_total = 0;
    int fb_total_t = 0;
    std::vector<Operation> fb_total_history;
};

// Experimental objective (off unless a JoinTable is passed). A t-BFB cycle
// appending a source that starts with segment v to a sequence that ends
// with segment u creates the join u -> v; is_foldback(u, v) says whether
// that join is fold-back-shaped. With a table, the distance is minimum
// first in the number of t-BFB cycles creating such a join, then in t-BFB
// cycles, then in the total number of cycles.
struct JoinTable {
    int max_id = 0;                      // identifiers lie in -max_id..max_id
    std::vector<std::uint8_t> foldback;  // (2 max_id + 1)^2 flags

    bool is_foldback(int u, int v) const {
        const std::size_t width = 2 * static_cast<std::size_t>(max_id) + 1;
        return foldback[static_cast<std::size_t>(u + max_id) * width +
                        static_cast<std::size_t>(v + max_id)] != 0;
    }
};

class MemoryLimitExceeded : public std::runtime_error {
  public:
    explicit MemoryLimitExceeded(const std::string& message)
        : std::runtime_error(message) {}
};

// Byte budget shared by every search of a process. Each search accounts
// for the tables it allocates and releases them when it is destroyed, so
// the budget bounds the memory of all worker threads together.
class MemoryBudget {
  public:
    explicit MemoryBudget(std::size_t limit_bytes);

    // Throws MemoryLimitExceeded when the budget would be exceeded.
    void acquire(std::size_t bytes);
    void release(std::size_t bytes);

    std::size_t limit() const { return limit_; }
    std::size_t used() const { return used_.load(); }
    std::size_t peak() const { return peak_.load(); }

  private:
    std::size_t limit_;
    std::atomic<std::size_t> used_;
    std::atomic<std::size_t> peak_;
};

// The sources and their reverse complements. The only data that depends
// on the sources alone, so one instance is shared by every target of a
// batch run.
struct PreparedSources {
    std::vector<std::vector<int>> forward;
    std::vector<std::vector<int>> reverse;

    PreparedSources() = default;
    explicit PreparedSources(std::vector<std::vector<int>> sources);
};

// Computes the nested BFB distance of target over the given sources.
// Throws std::invalid_argument when a source is empty (as the Python
// validate_sources does) or when the instance exceeds the index limits,
// and MemoryLimitExceeded when the budget is exhausted.
Result nbfb_fast_dist(const PreparedSources& sources,
                      const std::vector<int>& target,
                      bool want_history,
                      MemoryBudget& budget,
                      const JoinTable* foldback_first = nullptr);

// Convenience overload with an unbounded budget.
Result nbfb_fast_dist(const std::vector<std::vector<int>>& sources,
                      const std::vector<int>& target,
                      bool want_history = false);

}  // namespace nbfb

#endif  // NBFB_FAST_HPP
