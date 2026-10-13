# `nbfb`: nested BFB distance in C++17

Computes the NBFB distance of a target over a source set, where the source set
is all contiguous substrings, in both orientations, of the given source lines.
The search runs over index tuples `E(phi) = t[1,phi]` and
`W(phi,a,b) = t[1,phi] D[a,b]`, with `D = rc(t) # g_1 # ... # g_r # rc(g_r) # ... # rc(g_1)`.
The successors of a node come in O(n) intervals, one occurrence query per cut
point, for O(n^4) time in the worst case. `PREFIX_FAST_DESIGN.md` gives the rules.
The answer is minimum first in t-BFB cycles and then in total cycles, or
fold-back-first with `--foldback-first`.

No external dependencies. The library is `nbfb.hpp/.cpp` (shared types),
`nbfb_blocks.hpp/.cpp` (the search) and `nbfb_state.hpp/.cpp` (the reference
string search). The command-line tool is `nbfb_cli.cpp`.

## Build

```
make
```

or

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
```

## Usage

### One instance

```
nbfb -i INSTANCE.txt [--history] [--max-mem-gb G]
```

`INSTANCE.txt` has one source per line and the target on the last line. Segments
are joined by `;` and written `chr0:1000-1900[R]` (forward) or `chr0:1900-1000[C]`
(reverse complement). The source set is the closure of the source lines under
contiguous substrings and reverse complements. A file that lists every block
explicitly therefore denotes the same set as a file that lists only the full
strings.

The last line of stdout is the total distance, an integer or `inf`. With
`--history`, three lines precede it:

```
t_cycles: 3
types: ttiiiiiiiii
operations: t,2;t,84;i,19,21;...
11
```

`types` is the sequence of cycle kinds of one optimal history. `operations` is
the same history: `p,x,y`, `i,x,y` or `f,x,y` with `P = v[:x]`, `Q = v[x:y]`,
`R = v[y:]`, or `t,k` with `k` the index of the appended block among the source
lines (`-1` if the block is not listed verbatim).

### Batch

```
nbfb --sources SOURCES.txt --targets TARGETS.txt [-j N] [-o OUT.tsv] [--history] [--max-mem-gb G]
```

One source set shared by many targets, one target per line. The tool writes one
TSV row per target, with the columns `index`, `tokens`, `t_cycles`, `total`,
`types`, `seconds`, and `ops` (with `--history`). `-j N` runs `N` worker threads.
Rows are written in input order. A target that fails (for example on the memory
budget) gives an `error` row, and the exit status is 1.

### Fold-back-first objective

```
nbfb -i INSTANCE.txt --foldback-first [--foldback-span BP|any]     # either mode
```

A t-BFB cycle appending a block that starts with segment `v` to a sequence
ending with segment `u` creates the join `u -> v`. The join is fold-back-shaped
when both segments lie on the same chromosome in opposite orientations, with
breakpoints at most `--foldback-span` bp apart (default 50,000; `any` drops the
distance test). The distance is then minimum first in such t-BFB cycles, then
in t-BFB cycles, then in total cycles. The flag adds a `foldback_t` line or
column, and such cycles print as `T`.

### Reference string search

```
nbfb -i INSTANCE.txt --string-search [--foldback-first] [--history]   # either mode
```

This is a direct search over strings, processed in order of length and keyed by
the string, with every cycle enumerated (`nbfb_state.cpp`). It takes the sources
exactly as listed, with no closure under substrings, so a comparison needs an
instance that lists every block. It is the exact reference for testing and is
practical up to about a hundred segments.

## Memory

All tables go through a shared byte budget. `--max-mem-gb G` (default 2) sets
the limit for the whole process. A search that would exceed it stops with exit
status 3. The largest tables are:
* the LCE table over the text, `4 (N+2)^2` bytes with `N = 2n + 2K + r + 1`
  (`K` the total length of the strings);
* the memoized occurrence queries;
* the union-find families, one per reached `(mask, phi, start)`.

## Verification

* **Replay.** Every returned history is replayed on strings inside the search.
  Each cycle must produce exactly the string of its tuple, every tuple must
  satisfy its invariants, and every appended block must be a source.
* **Agreement with the string search.** I compared the two searches on 35 real
  instances (the CASTLE loci and the COLO829 inputs), under both objectives,
  with the string search run over the explicit closure. All 70 results were
  identical.

Timings, one core:

| instance | n | sources | `nbfb` | `--string-search` |
|---|---:|---:|---:|---:|
| HCC1395 chr6q minimal | 81 | 552 | 0.3 s | 63 s |
| HCC1395 chr6q | 71 | 552 | 0.1 s | not run |
| HCC1954 chr15q | 87 | 362 | 1.3 s | not run |
