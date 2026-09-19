# `nbfb` — C++17 port of `NBFB_fast_dist`

An exact port of `src/bfb_dist/nbfb_fast.py` (`NBFB_fast_dist`): the breadth-first
search over index tuples of the paper appendix "An O(2^m n^4)-Time Algorithm for
NBFB Distance", with the state extended by the number of t-BFB cycles applied so
far, so that the answer is minimum first in t-BFB cycles and then in the total
number of cycles. Sources may be reused, so the state is the pair (t-BFB count,
index tuple) rather than a used-set.

No external dependencies. Two files make the library (`nbfb_fast.hpp` /
`nbfb_fast.cpp`) and one makes the command-line tool (`nbfb_cli.cpp`).

## Build

```
make                      # g++ 15.2.0, -std=c++17 -O2 -Wall -Wextra, -pthread
```

or

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build
```

Both produce the `nbfb` binary (the Makefile in this directory, CMake in
`build/`). Tested with g++ 15.2.0 and cmake 4.2.3.

## Usage

### One instance

```
nbfb -i INSTANCE.txt [--history] [--max-mem-gb G]
```

`INSTANCE.txt` is the format `bfb_dist.__main__.read_instance` reads: one source
per line, the target on the last line, segments joined by `;` and written
`chr0:1000-1900[R]` (forward) or `chr0:1900-1000[C]` (reverse complement). The
last line of stdout is the total distance, an integer or `inf`, exactly as
`python3 -m bfb_dist distance --mode nbfb-fast` prints it. With `--history`,
three lines precede it:

```
t_cycles: 5
types: ttfffiifiifitiifffifiifiiitt
operations: t,412;t,55;f,0,1;f,1,2;...
28
```

`types` is the sequence of cycle kinds of one optimal history, `operations` the
same history in the serialization of `bfb_dist.history.Operation` (`p,x,y`,
`i,x,y`, `f,x,y` or `t,source`), semicolon separated, replayable with
`State.apply`.

### Batch

```
nbfb --sources SOURCES.txt --targets TARGETS.txt [-j N] [-o OUT.tsv] [--history] [--max-mem-gb G]
```

One source set (one source per line) shared by many targets (one target per
line, same token format). One TSV row per target:

| column | meaning |
|---|---|
| `index` | 0-based position in `TARGETS.txt` |
| `tokens` | number of segments in the target |
| `t_cycles` | number of t-BFB cycles, or `inf`, or `error` |
| `total` | total number of cycles, or `inf`, or `error` |
| `types` | cycle kinds of one optimal history, e.g. `tiftf` (the error message when the row failed) |
| `seconds` | wall time for that target |
| `ops` | the full history, only with `--history` |

`-j N` runs `N` `std::thread` workers over the targets (default 1); rows are
always written in input order. A target that exhausts the memory budget or
overflows the index limits is reported as an `error` row and the run continues;
the exit status is then 1.

### Library

```cpp
#include "nbfb_fast.hpp"

nbfb::MemoryBudget budget(2ull << 30);
nbfb::PreparedSources sources({{1, 2}, {-2, -1}});   // signed segment ids
nbfb::Result r = nbfb::nbfb_fast_dist(sources, {1, 2, -2, -1}, true, budget);
// r.reachable, r.t_cycles, r.total, r.history
```

A sequence is a vector of signed segment identifiers: every distinct segment
interval gets a positive identifier and its reverse complement the negated one,
so reverse complementing a sequence is reversing the order and negating every
entry. The CLI builds that map from the token strings (two tokens get the same
magnitude exactly when one is the `Segment.reverse_complement()` of the other,
matching the Python `Segment` equality on `(chr, start, end, complement)`).

## Design notes

### Index layout

All target indices are 1-based to mirror the paper, and every table is a flat
`std::vector` indexed arithmetically with a 32-bit element type.

* **LCE tables.** `LCE(t,t)`, `LCE(tbar,t)` and `LCE(t,tbar)` are stored in
  full, `(n+2) x (n+2)` `int32` each, row-major with stride `n+2`; the padding
  row and column let a query run one position past either end, as in the Python.
* **Source LCE tables.** The algorithm only ever reads row 1 of `LCE(x,t)` and
  column 1 of `LCE(t,xbar)`, so only those slices are stored, `m x (n+2)`
  `int32` each (the Python builds both tables in full). The occurrence table
  does need all of `LCE(t,x)`; that one is built lazily, only for the sources
  that actually appear as a source node `T(phi0, x)`.
* **OCC / MINB.** Memoized in two open-addressing hash maps (linear probing,
  power-of-two capacity, `splitmix64`), keyed by the packed 64-bit
  `zkey<<48 | i<<32 | lambda<<16 | s` with the value pair `(best, pos)` as two
  `int32`. The dense index space is `O(m n^3)`, while an instance only queries a
  small part of it, so a hash map is the right structure here; the search itself
  uses arrays instead (below). The maps are built once per instance and survive
  the deepening loop over the t-BFB bound, exactly as the Python dictionaries do.
* **Visited set.** Per layer, that is per value of the t-BFB count:
  * `E(phi)`: one interval union-find over `0..n`, `(n+2)` `int32`.
  * `W(phi,a,b)`: one interval union-find over `a in 1..n` for each `(phi,b)`
    pair. A per-layer index array of `(n+1)^2` `int32` maps `(phi,b)` to a block
    id, and blocks of `n+1` `int32` are appended to a pool the first time a
    `(phi,b)` pair is reached.
  * `T(phi0,x)`: one byte per `(phi0, x)`, `(n+1) * m` bytes, allocated on first
    use.
* **Queue, distances, parents.** One flat `std::vector<Record>` of 32-byte
  records holds every visited tuple in discovery order. The breadth-first queue
  is a head index into that vector (records are never removed), the distance is
  a field of the record, and the parent pointer plus the operation that produced
  the record reconstruct one optimal history. There is no per-node heap
  allocation and no hash map anywhere in the search loop.

### Memory bound and `--max-mem-gb`

Every allocation above goes through a shared `MemoryBudget`, an atomic byte
counter. `--max-mem-gb G` sets the limit for the whole process (all `-j`
workers together); an allocation that would exceed it throws, and the tool
prints `Error: memory limit exceeded ... (raise --max-mem-gb to continue)` and
exits with status 3 rather than growing. The default is 2 GB.

Accounted bytes, with `n` the target length and `m` the number of sources:

| table | bytes |
|---|---|
| three full LCE tables | `12 (n+2)^2` |
| source LCE slices | `8 m (n+2)` |
| lazy `LCE(t,x)` | `4 (n+2)(len(x)+2)` per source actually used |
| `keep_w` / `keep_t` | 12 per entry, at most `12 n m` |
| E family, per layer | `4 (n+2)` |
| W index, per layer | `4 (n+1)^2` |
| W blocks, per layer | `4 (n+1)` per reached `(phi,b)` pair, at most `~4 n^3` |
| T flags, per layer | `(n+1) m` |
| records | 32 per visited tuple |

The number of layers is one more than the minimum number of t-BFB cycles (the
deepening loop rebuilds the layers for each bound, so only the last bound's
layers are live).

The W family is the term the task singled out: a dense `(n+1) x n x (n+1)` array
is `1.6e7 int32 = 64 MB` per layer at `n = 250` even when the search touches a
handful of `(phi,b)` pairs, while a hash map would cost a hash per `next` query
in the innermost loop of the interval union-find. The two-level layout pays
`4 (n+1)^2` bytes up front per layer (250 KB at `n = 250`), keeps one indirection
per union-find access, and reaches the `4 n^3` worst case (62 MB per layer at
`n = 250`) only if every `(phi,b)` pair is actually visited.

Measured peak RSS (`/usr/bin/time -v`, whole target file, `--max-mem-gb 3`):

| locus | tokens | sources | `-j 1` | `-j 4` |
|---|---:|---:|---:|---:|
| HCC1954 chr11q (2,000 walks) | 43 | 56 | 6.9 MB | 7.4 MB |
| HCC1395 chr6q (240 walks) | 89 | 552 | 42 MB | 101 MB |
| HCC1954 chr15q P.cen (240 walks) | 81 | 364 | 23 MB | 63 MB |
| H1437 chr5p Q.cen (60 walks) | 203 | 514 | **206 MB** | 794 MB |

The budget counts the search tables, not allocator overhead or the CLI input
buffers, so the resident set sits somewhat above the accounted total; set
`--max-mem-gb` below the real ceiling of the machine.

### What is shared across targets in batch mode

Only the data that depends on the sources alone: the segment-identifier map and
the `PreparedSources` (the sources and their reverse complements). Every LCE,
occurrence and `keep_*` table involves the target, so those are rebuilt per
target, and their memory is returned to the budget when the target is done.

## Verification

`verify.py` is the differential harness. It solves every instance twice, once
with `bfb_dist.nbfb_fast.NBFB_fast_dist` and once with `nbfb`, requires the total
distance *and* the t-BFB count to agree, and replays the C++ history with the
Python `State.apply` to check that it reaches the target.

```
make
PYTHONPATH=../src python3 verify.py --part all -n 3000 --unreachable-count 1000 \
    --reference-count 500 -j 4 --time-python 60 --max-mem-gb 3 \
    --scratch /path/to/scratch
```

Parts: `random` (simulated histories), `unreachable` (random targets),
`reference` (the tiny instances are also solved with the `O(n^6)`
`bfb_dist.nbfb.NBFB_dist`, so a wrong successor rule shared by the Python and
the port would show up), `castle` (the stored longest walks), `all`.

### Results (2026-09-19, g++ 15.2.0, `-O2`)

| instances | checked | mismatches | history replays |
|---|---:|---:|---|
| simulated histories (`simulate_history`, 1-4 sources, 3-12 segments, 1-8 cycles, half with `short_deletions`, one third with the source set replaced by all contiguous blocks of the sources in both orientations) | 3,000 | **0** | 2,872 / 2,872 (the other 128 are unreachable) |
| random targets over the source alphabet, some with a segment absent from every source | 1,000 | **0** | 103 / 103 (897 of the 1,000 unreachable) |
| tiny instances also solved with the `O(n^6)` `NBFB_dist` (C++ = fast Python = reference) | 500 | **0** (0 fast-vs-reference) | 499 / 499 |
| HCC1954 chr11q longest walks | 2,000 | **0** | 2,000 / 2,000 |
| HCC1954 chr11q longest walks (extra run) | 20,000 | **0** | 20,000 / 20,000 |
| HCC1395 chr6q longest walks | 240 | **0** | 240 / 240 |
| HCC1954 chr15q P.cen longest walks | 240 | **0** | 240 / 240 |
| H1437 chr5p Q.cen longest walks | 60 | **0** | 60 / 60 |

The CASTLE rows are checked against the `nbfb` and `nbfb_types` fields stored in
`CASTLE/derived/longest_walks/*.json` (total distance and, from `nbfb_types`, the
number of t-BFB cycles). The instance is rebuilt exactly as
`CASTLE/scripts/longest_walks.py` builds it (`_init` / `_tok`: chains from
`molecules.build_chains` with the model and `--append` arguments of
`CASTLE/derived/LONGEST_WALKS.md`, synthetic coordinates, sources = every
contiguous block of every chain in both orientations).

On the walks tried, the C++ reports not only the same distance but the same
cycle-type string as the stored Python history, because the tie-breaking of the
successor order was kept identical (`std::stable_sort` by decreasing extension
length over the sources in increasing index order, first occurrence of `xbar`).

The existing Python unit tests still pass:

```
PYTHONPATH=src python3 -m unittest -q tests.test_nbfb tests.test_source_use tests.test_cli
# Ran 29 tests ... OK
```

The two command-line tools were also compared directly
(`python3 -m bfb_dist distance --mode nbfb-fast` versus `nbfb -i`), and batch
mode against single-instance mode, on chr11q instances: identical last lines,
identical t-BFB counts and identical cycle types.

### Speed, single core

`--time-python` re-times Python and the C++ on the same walks in the same run,
one worker each.

| locus | tokens | sources | Python | C++ | speedup |
|---|---:|---:|---:|---:|---:|
| chr11q | 43 | 56 | 5.4 ms/walk | 0.1 ms/walk | 58x |
| chr6q | 89 | 552 | 518 ms/walk | 10.9 ms/walk | 48x |
| chr15q P.cen | 81 | 364 | 636 ms/walk | 18.4 ms/walk | 35x |
| chr5p Q.cen | 203 | 514 | 12.1 s/walk | 396 ms/walk | 31x |

Whole target file, wall clock:

| locus | walks | `-j 1` | `-j 4` |
|---|---:|---:|---:|
| chr11q | 2,000 | 0.19 s | 0.06 s |
| chr6q | 240 | 2.04 s | 0.58 s |
| chr15q P.cen | 240 | 4.43 s | 1.19 s |
| chr5p Q.cen | 60 | 23.3 s | 5.96 s |

## Notes on the Python

Nothing that changes an answer was found: the port is a faithful translation and
every instance agrees. Three observations about `src/bfb_dist/nbfb_fast.py`:

1. `self.z_len` is filled in `_FastSearch.__init__` and never read. `self.ext_w`
   and `self.ext_t` are read only while building `keep_w` and `keep_t`, so they
   could be local as well; on a CASTLE-sized instance they are `m` rows of
   `n + 2` tuples that stay alive for the whole search.
2. `self.lce_tx` and `self.lce_txbar` are built as full `(n+2) x (|x|+2)` tables
   for every source, but only row 1 of `lce_xt` and column 1 of `lce_txbar` are
   ever read (`lce_xt[index][1][phi+1]`, `lce_txbar[index][p][1]`). On chr5p
   (514 sources) this is the dominant table cost. The port stores only the two
   slices and builds `LCE(t,x)` lazily.
3. In `expand_tail`, the expression `self.keep_w[a] if node[0] == "W" else
   self.keep_t[index]` relies on Python evaluating only the selected branch: `a`
   is undefined on the source-node path and `index` is undefined on the
   window-node path. Rewriting it as a plain `if`/`else` with both values
   computed would raise `NameError`.

One structural remark, shared by both implementations: the index-tuple search has
no analogue of the `len(child) > n` pruning of `nbfb._successors`. A window node
`W(phi,a,b)` represents a string of length `phi + (b - a + 1)`, which can exceed
`n`. Such a node can never reach `E(n)` (the guard `lo <= hi` on the empty-tail
interval and the monotone length growth of the successor rules make sure of it),
so the distance is unaffected, but the search does expand tuples representing
strings longer than the target.
