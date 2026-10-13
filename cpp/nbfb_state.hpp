// Nested BFB distance by a direct search over strings: the exact
// reference for nbfb_blocks_dist.
//
// Every BFB cycle lengthens the string, so the graph whose nodes are the
// strings of length at most |t| and whose edges are single p-, i-, f- and
// t-BFB cycles is acyclic, and a minimum-cost path is found by processing
// the strings in order of length. Every produced string must be nested:
// t[1,psi] y with y empty, the reverse complement of a substring of t, or a
// source. The search keys nodes by the string itself and enumerates every
// cycle explicitly, so it is O(n^2) cycles per node and is meant for
// testing on instances of up to a few hundred segments.

#ifndef NBFB_STATE_HPP
#define NBFB_STATE_HPP

#include <vector>

#include "nbfb.hpp"

namespace nbfb {

// Minimum first in fold-back-shaped t-BFB joins (only with a JoinTable),
// then in t-BFB cycles, then in total cycles. The sources are taken as
// given (no closure under substrings).
Result nbfb_state_dist(const PreparedSources& sources,
                       const std::vector<int>& target,
                       bool want_history,
                       MemoryBudget& budget,
                       const JoinTable* foldback_first = nullptr);

}  // namespace nbfb

#endif  // NBFB_STATE_HPP
