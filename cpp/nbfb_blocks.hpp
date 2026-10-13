// Nested BFB distance when the sources are all contiguous substrings, in
// both orientations, of a few strings g_1, ..., g_r.
//
// The given source lines are read as those strings: the source set is the
// set of all their substrings and the reverse complements of those. A file
// that lists every block explicitly therefore denotes the same set.
//
// Every node is t[1,phi] D[a,b] with D = rc(t) # g_1 # ... # g_r #
// rc(g_r) # ... # rc(g_1), and the successors of a node come in O(n)
// intervals of index tuples. See PREFIX_FAST_DESIGN.md.

#ifndef NBFB_BLOCKS_HPP
#define NBFB_BLOCKS_HPP

#include <vector>

#include "nbfb.hpp"

namespace nbfb {

// Minimum first in fold-back-shaped t-BFB joins (only with a JoinTable),
// then in t-BFB cycles, then in total cycles. For a t-BFB cycle,
// Operation::source is the index of the appended block among the given
// sources, or -1 if the block is not listed verbatim.
Result nbfb_blocks_dist(const PreparedSources& sources,
                        const std::vector<int>& target,
                        bool want_history,
                        MemoryBudget& budget,
                        const JoinTable* foldback_first = nullptr);

}  // namespace nbfb

#endif  // NBFB_BLOCKS_HPP
