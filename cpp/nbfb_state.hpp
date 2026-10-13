// Nested BFB distance by a direct search over strings, for comparing
// definitions of nestedness.
//
// Every BFB cycle lengthens the string, so the graph whose nodes are the
// strings of length at most |t| and whose edges are single p-, i-, f- and
// t-BFB cycles is acyclic, and a minimum-cost path is found by processing
// the strings in order of length. The search keys nodes by the string
// itself and enumerates every cycle explicitly, so it is O(n^2) cycles per
// node and is meant for instances of up to a few hundred segments, not as
// a replacement for nbfb_fast_dist.
//
// Two definitions restrict the nodes and edges:
//
//   Definition::Prefix   (nested since 2026-10-02) every produced string is
//       t[1,psi] y with y empty, y the reverse complement of a substring of
//       t, or y a source. No condition on the cycle.
//   Definition::Lcp      (nested until 2026-10-02) every p-, i- and f-BFB
//       cycle strictly extends the longest common prefix with t, and every
//       produced string has a foldable tail unless it is produced by a
//       t-BFB cycle that extends the longest common prefix. This is the
//       class searched by nbfb_fast_dist and is kept here for differential
//       testing.

#ifndef NBFB_STATE_HPP
#define NBFB_STATE_HPP

#include <vector>

#include "nbfb_fast.hpp"

namespace nbfb {

enum class Definition { Prefix, Lcp };

// Same result conventions as nbfb_fast_dist: minimum first in fold-back-
// shaped t-BFB joins (only with a JoinTable), then in t-BFB cycles, then in
// total cycles. The fb_total fields are not filled.
Result nbfb_state_dist(const PreparedSources& sources,
                       const std::vector<int>& target,
                       bool want_history,
                       MemoryBudget& budget,
                       Definition definition,
                       const JoinTable* foldback_first = nullptr);

}  // namespace nbfb

#endif  // NBFB_STATE_HPP
