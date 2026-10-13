# Fast NBFB search for the prefix definition of nestedness (rev. 4, implemented in nbfb_blocks.cpp)

Goal: an index-tuple search for `def:nested` with a provable O(n^4) bound.
In `def:nested`, every string is `t[1,psi] y`, with `y` empty, the reverse
complement of a substring of `t`, or a source. `nbfb --string-search`
(the string search) stays as the exact oracle for testing. The LCP-monotone
class is dropped.

Scope (author decision, 2026-10-13): given strings g_1, ..., g_r (one per
reference-contiguous run of segments, as in every real-data analysis),

    S = all contiguous substrings of the g_i and their reverse complements.

Rev. 4 replaces `t` by the g_i in rev. 3: `D = tbar # g_1 # ... # g_r #
rc(g_r) # ... # rc(g_1)` (the code appends `# t` so that the text is its own
reverse complement, and restricts occurrences to D). A covering source ends
with the tail, so the tail is a substring of some g_i or rc(g_i). Case C
appends a separator-free window of the g part. `nbfb` reads the
source lines as the g_i (a file listing every block denotes the same S).

Notation:
* `n = |t|`.
* `rc(.)` is the reverse complement, and `tbar = rc(t)`.
* `phi` is the length of the longest common prefix of a string with `t`.
* `[x,y]` is an inclusive range.

## 1. Nodes: a single tuple type

Let `D = tbar # t`, with `N = |D| = 2n+1`. The separator `#` compares
unequal to every symbol, including itself.

**Lemma (node test).** A string `t[1,phi] u`, with `phi` exact, is a node
iff `u` occurs in `D`. Since `u` contains no `#`, this means `u` occurs in
`tbar` or in `t`.

*Proof.*
* By the proof of `thm:nbfb`, the string is a foldable prefix iff `rc(u)`
  is a substring of `t`, that is, iff `u` occurs in `tbar`.
* By the same proof, the string is `t[1,psi] y` with `y` in `S` iff some
  source `y` is a suffix of the string with `|y| >= |u|`.
  * If such a `y` exists, `u` is a suffix of `y`. So `u` is a substring of
    `t`, which means `u` occurs in `t`.
  * Conversely, if `u` occurs in `t`, then `u` is a source. Take `y = u` and
    `psi = phi`.
* Every source occurs in `t`, so a t-BFB appends exactly the strings
  `t[i,j]`.

The tuples are:
* `E(phi)`: `t[1,phi]`.
* `W(phi, a, b)`: `t[1,phi] D[a,b]`, with `D[a,b]` free of `#` and `phi`
  exact, that is `D[a] != t[phi+1]`.

There are O(n N^2) = O(n^3) tuples. Several tuples can denote the same
string, which is harmless, as in the old appendix. Every tuple has a known
length, so children longer than `n` are clipped, including inside
intervals.

## 2. One formula for every fold-back child

Let `s` be the current string. Every p-, i-, and f-BFB child is

    child(c, beta, gamma) = s[1,c] . rc(s[beta+1, gamma])

* p/i cycles: `c = |s|` and `0 <= beta < gamma <= |s|`. The cycle is a
  p-BFB when `gamma = |s|`.
* f cycles: `gamma = |s|` and `0 <= beta < c < |s|`.

## 3. The monotone structure behind every rule

Fix the cut point (`c`, or `gamma`) and the matched length. Every case
below then has the same shape.

* The child's tail is `tail(beta) = L . G(beta)`.
* `L` is a fixed left part.
* `G(beta)` grows at its right end as `beta` decreases.
* Hence `tail(beta)` occurs in `D` for every `beta` in an interval
  `[beta_min, beta_max]`.

**One query answers each cut point.**

* `RightExt(L, G)` is the maximum, over occurrences `p` of `L` in `D`, of
  `LCE(D[p+|L|, N], G)`.
* `beta_min` follows from `RightExt`.
* Let `p*` be an argmax occurrence. Every valid `tail(beta)` is then a
  prefix of `D[p*, ...]`, so the children are
  `W(phi', p*, p* + |tail(beta)| - 1)`.
* These children form one interval in the union-find family `(phi', p*)`,
  whose free index is the end `b'`.

**Computing `RightExt` in O(1).**

* Mirror the problem: `tail` occurs in `D` iff `rc(tail)` occurs in
  `U = rc(D)`. Since `rc(D) = tbar # t`, we have `U = D` up to separator
  placement. An occurrence ending at `q` in `U` starts at `N - q + 1` in
  `D`.
* `RightExt` becomes a left extension of the fixed right part `rc(L)`. It
  is a range-max, over the suffix-array interval of `rc(L)`, of the longest
  common suffix of `U[1, p-1]` with the context. That value does not depend
  on `|L|`, so one range-max table per context column suffices. The tables
  store argmax positions as well as values.
* Contexts end in a prefix of `t` or in a window of `D`. They need LCS
  tables `U x t` and `U x D`, (rank, length) suffix-array intervals over a
  generalized suffix array of `U` and `D`, and sparse tables. That is
  O(n^2 log n) time and memory.
* When `G` consists of two pieces, the query runs in two stages:
  1. Match the first piece.
  2. Only if the first piece is matched completely, take the suffix-array
     interval of the longer fixed part and match the second piece.

## 4. Successor rules

Let `s = t[1,phi] u` with `u = D[a,b]` and `lambda = |u|`.

| case | cycles | `L` (fixed) | `G(beta)` (grows) |
|---|---|---|---|
| E node | all p/i/f; t-BFB | (direct) | every child is a node. Fold-backs are enumerated directly; t-BFB children are `W(phi, t-part of D: i, j)`, as one interval in `j` per start `i` |
| A full | f, `1 <= c <= phi`, `rc(u)` matched completely by `t[c+1..]` | (direct) | old full-absorption intervals with the guards dropped |
| A partial | f, `1 <= c <= phi`, `M < lambda` symbols of `rc(u)` matched | `rc(u)[M+1, lambda]` | `rc(t[beta+1, phi])` |
| B1 | f, `c = phi+k`, `1 <= k < lambda` | `u[1,k]` | `rc(u)` then `rc(t[beta+1, phi])` (two stages) |
| B2 | i, `gamma <= phi` | `u` | `rc(t[beta+1, gamma])` |
| B3 | p/i, `gamma = phi+g`, `1 <= g <= lambda` | `u` | `rc(u[1,g])` then `rc(t[beta+1, phi])` (two stages) |
| C | t-BFB appending `t[i,j]` | `u` | `t[i,j]`, which grows in `j` for fixed `i` |

**Case A.**
* The full-absorption intervals are the old ones without the "`phi` must
  grow" guards.
* Partial absorption runs for every `M` in `[0, lambda-1]`, through
  `RightExt` (no OCC/MINB tables).

**Common prefix of the children.**
* In cases B and C, `phi` is unchanged, because `u[1]` still mismatches
  `t[phi+1]`.
* In case A, the child's common prefix is `c+M`, which is exact.
* Only E nodes and full absorption produce `E` children.

**Case C.** For each start `i` in `[1,n]`, one query gives the largest `j`.
That is O(n) queries per node. The extension stops at `#` and at the end of
`D`, so every appended string is a substring of `t`.

The case table was reviewed for rev. 2 (general substring-closed `S`) and
found complete and sound. Rev. 3 is the special case `Z = t`.

## 5. Search and complexity

* Search: the appendix's single lexicographic bucket-Dijkstra over tuples,
  with union-find families keyed by `(phi', a')` and free index `b'`. No
  iterative deepening is needed.
* Expanding a reachable `W` node costs O(|s| + n) = O(n) queries plus O(1)
  amortized per visited child.
* Tables cost O(n^2 log n).
* Total: O(n^2 log n + sum over reachable tuples of n), which is O(n^4) in
  the worst case and output-sensitive in the reachable tuples.

## 6. Verification

* Replay every returned history on strings.
* Run differential tests against the string search, with `S` = all
  substrings of `t`.
* Disable each rule in turn and confirm that the oracle catches it.

## 7. Input

* `nbfb --substrings -i FILE` reads only the target (the last line) and uses
  `S` = all substrings of `t`.
* Duplicate tuples for the same string are allowed, as in the old appendix.
