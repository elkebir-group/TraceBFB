# `nbfb`: nested BFB distance in C++17

Computes the NBFB distance of a target over a source set, one optimal history,
and the number of optimal histories, and optionally lists every optimal history. A history is nested if every string it produces is a foldable
prefix of the target `t`: `s = x y`, where `x` is a prefix of `t` and `y` or
`rc(y)` is a substring of `t`.

No external dependencies. The library is `nbfb.hpp/.cpp` (shared types) and
`nbfb_simple.hpp/.cpp` (the search). The command-line tool is `nbfb_cli.cpp`.

## Build

```
make
```

or

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
```

## Usage

```
nbfb -i INSTANCE.txt [options]                                     # one instance
nbfb --sources S.txt --targets T.txt [-j N] [-o OUT.tsv] [options] # batch
```

`INSTANCE.txt` has one source per line and the target on the last line. Segments
are joined by `;` and written `chr0:1000-1900[R]` (forward) or `chr0:1900-1000[C]`
(reverse complement). In batch mode every line of `T.txt` is a target over the
shared sources of `S.txt`.

| option | meaning |
|---|---|
| `--blocks` (default) | sources = every contiguous block of the source lines and its reverse complement (a file listing every block denotes the same set) |
| `--general` | sources = the source lines exactly as listed |
| `--foldback-first` (default) | minimize fold-back-shaped t-BFB joins, then t-BFB cycles, then all cycles |
| `--no-foldback-first` | minimize t-BFB cycles, then all cycles |
| `--foldback-span BP\|any` | a join `u -> v` is fold-back-shaped when both segments lie on the same chromosome in opposite orientations with breakpoints at most `BP` apart (default 50,000; `any` drops the distance test) |
| `--all-optimal` | list every optimal history |
| `--max-histories N` | list at most `N` optimal histories |
| `--max-mem-gb G` | memory budget for the whole process (default 2); exceeding it stops with exit status 3 |
| `-j N`, `-o FILE` | batch mode: worker threads and output TSV (default stdout) |

### Output

Single mode: the last line of stdout is the distance, an integer or `inf`. Before
it come `foldback_t:` (with `--foldback-first`), `t_cycles:`, `types:`,
`operations:`, `optimal_histories: COUNT`, and one `optimal: TYPES<tab>OPERATIONS`
line per listed history (with `--all-optimal` or `--max-histories`).

Batch mode: one TSV row per target with the columns `index`, `tokens`,
`foldback_t` (with `--foldback-first`), `t_cycles`, `total`, `types`, `seconds`,
`ops`, `optimal_histories`, and `all_ops` (with `--all-optimal` or
`--max-histories`; histories separated by `|`). A target that fails, for
example on the memory budget, gives an `error` row, and the exit status is 1.

`types` lists the cycle kinds, with `T` for a fold-back-shaped t-BFB cycle.
`operations` gives `p,x,y`, `i,x,y` or `f,x,y` with `P = v[:x]`, `Q = v[x:y]`,
`R = v[y:]`, or `t,k` with `k` the index of the appended source among the source
lines (`-1` if `--blocks` appends a block that is not listed verbatim).

Two histories are the same when they have the same strings and the same cycle
types. For each step of a listed history, one cut point or source is printed.

## Algorithm

A Dijkstra search over the foldable prefixes of `t`, each stored as
`t[1,phi] T[c, c+l)` with `T = t # rc(t)`. For a fixed cut point, the p-, i- and
f-children of a string are the prefixes of one of two strings in a range of
lengths. Foldable prefixes are closed under taking prefixes, so a binary search
finds the children that are nodes, and union-find pointers skip the children that
a cycle of the same cost class already reached. With `--general`, t-BFB children
come from a pruned walk over a trie of the sources. With `--blocks`, the maximal
source lines are appended to the text and the t-BFB children of one start
position form a chain like the fold children. Running time:
`O(n^4 log n + n^3 |S|)` with `--general` and `O(n^3 (n + K) log n)` with
`--blocks`, for `n = |t|` and `K` the total length of the maximal source lines.

The count comes from a backward pass from `t` that inverts every cycle and keeps
the steps `u -> v` with `d(u) + cost = d(v)`. Every string of an optimal history
has its exact cost, because every cycle raises the cost. The count is a path
count (128-bit), and the listing a depth-first search.

## Verification

* **Replay.** Every returned history, and every listed optimal history, is
  replayed on strings: each cycle must be valid, every string a foldable prefix
  of `t`, the last string `t`, and the cost the reported one.
* **Brute force.** `tests/test_cpp_simple.py` compares both source modes and both
  objectives with a Python search that applies the definition literally, including
  the exact set of optimal histories.

Timings, one core, `--foldback-first`:

| instance | n | sources | `--blocks` |
|---|---:|---:|---:|
| HCC1954 chr11q | 43 | 56 | 0.01 s |
| HCC1395 chr6q | 71 | 552 | 0.1 s |
| H1437 chr5p | 191 | 514 | 7–12 s |

# `lwalk`: longest walks of a breakpoint graph

```
lwalk GRAPH.txt [--starts LABEL,...] [--max-foldbacks] [--max-mem-gb G]          # count
lwalk GRAPH.txt ... --enumerate -o PREFIX [--dedup-rc] [--max-walks N] [--chunk N]
lwalk GRAPH.txt ... --sample N [--seed S] -o PREFIX [--dedup-rc]
```

A walk starts at a loose end of the graph and uses every segment at most its
copy number and every bond at most its multiplicity. `lwalk` counts the longest
walks exactly (128-bit), enumerates them, or samples them uniformly, by dynamic
programming over (start, current node, remaining bond multiplicities). With
`--max-foldbacks` it maximizes the number of fold-back junction traversals first
and the length second. When the graph is the graph of one molecule, the longest
walks are the walks that use every edge.

`PREFIX.walks.txt` gets one walk per line as segment names (lower case for a
reverse traversal) and `PREFIX.targets.txt` the same walks as `;`-joined tokens,
the target format of `nbfb --targets`. `tracebfb.breakpoint_graph.BreakpointGraph.write_lwalk`
writes `GRAPH.txt`.
