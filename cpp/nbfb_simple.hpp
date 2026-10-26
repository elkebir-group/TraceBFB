// Nested BFB distance over a general source set.
//
// A string s is a foldable prefix of the target t if s = x y, where x is a
// prefix of t and y or the reverse complement of y is a substring of t. A
// history is nested if every string of the history is a foldable prefix of
// t. The sources are taken exactly as listed (no closure under substrings).
//
// The search is a Dijkstra search over the foldable prefixes of t under the
// lexicographic cost. A node is stored as t[1,phi] T[c,c+l), where
// T = t # rc(t), phi is the longest common prefix of the node and t, and
// T[c,c+l) is the first occurrence of the tail in T.
//
// p-, i- and f-cycles: for a fixed cut point b, the children of a node w
// are the prefixes, in a range of lengths, of one string X (w rc(w[1,b])
// for p and i, w[1,b] rc(w[b+1,|w|]) rc(w[1,b]) for f). Foldable prefixes
// are closed under taking prefixes, so a binary search finds the children
// that are nodes. All these cycles cost the same, so the first time a fold
// reaches a node is the best fold offer, and union-find pointers over
// prefix lengths skip the children that a fold already reached. Two exact
// necessary conditions skip most chains before the binary search.
//
// t-cycles: a trie of the sources is walked from w and pruned at the
// first prefix that is not a foldable prefix.
//
// String equality is tested through hashing, and every hash match is
// confirmed by an exact comparison, so results do not depend on the hash.
// Every returned history is replayed on strings before it is reported.
// The running time is O(n^4 log n + n^3 |S|) for n = |t|, up to the trie
// walk; the memory is O(n^2) for the tables plus O(1) per reached node.

#ifndef NBFB_SIMPLE_HPP
#define NBFB_SIMPLE_HPP

#include <vector>

#include "nbfb.hpp"

namespace nbfb {

// Minimum first in fold-back-shaped t-BFB joins (only with a JoinTable),
// then in t-BFB cycles, then in total cycles. For a t-BFB cycle,
// Operation::source is the index of the appended source.
//
// With blocks, the source set is every contiguous block of the source lines
// g_1, ..., g_r and the reverse complement of every such block (a file that
// lists every block explicitly denotes the same set). For a t-BFB cycle,
// Operation::source is then the index of the appended block among the
// source lines, or -1 if the block is not listed verbatim. The t-BFB
// children for a fixed start position in some g_k or rc(g_k) are the
// prefixes of one string, so they form chains like the fold children, for
// O(n^3 (n + K) log n) time with K the total length of the maximal lines.
//
// The result gives one optimal history and the exact number of optimal
// histories, and lists up to max_histories of them (Result::optimal_count,
// Result::optimal_histories). Two histories are the same when they have the
// same strings and the same cycle types; one cut point or source stands for
// each step. The search keeps one parent per node, so a second, backward
// pass builds the DAG of all optimal steps: every string of an optimal
// history is settled with its exact cost, and a step u -> v is optimal when
// d(u) + cost = d(v).
Result nbfb_simple_dist(const std::vector<std::vector<int>>& sources,
                        const std::vector<int>& target,
                        MemoryBudget& budget,
                        const JoinTable* foldback_first = nullptr,
                        bool blocks = true,
                        std::size_t max_histories = 0);

}  // namespace nbfb

#endif  // NBFB_SIMPLE_HPP
