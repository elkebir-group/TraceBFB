# TraceBFB

TraceBFB reconstructs parsimonious breakage-fusion-bridge (BFB) histories with
deletions at fold-back junctions and translocations, from assembled tumor
sequences and from breakpoint graphs.

TraceBFB solves the problems of the paper *TraceBFB: Reconstructing
Parsimonious Deletion- and Translocation-Aware Breakage-Fusion-Bridge
Histories*:

| Problem | Input | TraceBFB |
|---|---|---|
| NBFB Distance | a set of sources and a target | `cpp/nbfb --general` |
| Interchromosomal NBFB Distance | a target | `cpp/nbfb` |
| NBFB Walk | a breakpoint graph and a start node | `cpp/lwalk`, then `cpp/nbfb --targets` |

It also provides:

- `tracebfb graph`, which builds the breakpoint graph of a locus from the
  long reads of a tumor sample and a matched normal sample.
- `tracebfb distance`, which computes the minimum cost of a history that need
  not be nested, by backward search (exponential time).
- `tracebfb simulate`, which simulates histories.

## Installation

The Python package requires Python 3.7 or newer and `pysam`:

```bash
python -m pip install .
```

This installs the `tracebfb` command and the `tracebfb` Python package.

The C++ tools `nbfb` and `lwalk` need a C++17 compiler and nothing else:

```bash
make -C cpp
```

With CMake instead:

```bash
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release && cmake --build cpp/build
```

## Model

Strings represent genomic fragments on one strand, and `rc(x)` is the reverse
complement of a string `x`. A set `S` of **sources** represents segments of the
reference genome, and a **target** `t` represents a segment of a rearranged
genome.

A BFB cycle transforms a string `s` into a longer string:

| Cycle | Name | `s` | Result | Condition |
|---|---|---|---|---|
| p-BFB | perfect BFB cycle | `p q` | `p q rc(q)` | `q` non-empty |
| i-BFB | BFB cycle with inverted-part deletion | `p q r` | `p q r rc(q)` | `q`, `r` non-empty |
| f-BFB | BFB cycle with forward-part deletion | `p q r` | `p q rc(r) rc(q)` | `q`, `r` non-empty |
| t-BFB | BFB cycle with translocation | `s` | `s x` | `x` a source in `S` |

The string `r` is the **overhang**. A **history** over `S` is a sequence of
p-, i-, f- and t-BFB cycles that transforms the empty string into the target
`t`. A source may be appended any number of times.

A string `s` is a **foldable prefix** of `t` if `s = x y`, where `x` is a prefix
of `t` and `y` or `rc(y)` is a substring of `t`. A history is **nested** if every
string of the history is a foldable prefix of `t`. Every history whose overhangs
are shorter than every duplicated string and every appended source is nested.

### NBFB Distance

Given a set `S` of non-empty sources and a target `t`, **NBFB Distance** seeks a
nested history over `S` that transforms the empty string into `t` with minimum
cost. `nbfb --general` appends the sources exactly as listed.

### Interchromosomal NBFB Distance

The **blocks** of `t` are the maximal contiguous parts of the normal genome that
the segments of `t` cover. The sources of a set of blocks are the non-empty
substrings of the blocks and the reverse complements of these substrings. A
t-BFB cycle `s -> s x` is a **fold-back t-BFB cycle** if `s` ends with a segment
`y` and `x` starts with `rc(y)`.

Given a target `t`, **Interchromosomal NBFB Distance** seeks a nested history
over the sources of the blocks of `t` that lexicographically minimizes the
cost `(ft, t, total)`. Here `ft` is the number of fold-back t-BFB cycles, `t` the
number of t-BFB cycles and `total` the number of all cycles.

`nbfb` solves Interchromosomal NBFB Distance by default. With `--blocks` (the
default), the source lines take the place of the blocks: a t-BFB cycle appends
any contiguous part of a source line or the reverse complement of that part.
With `--foldback-first` (the default), the cost is `(ft, t, total)`. With
`--no-foldback-first`, the cost is `(t, total)`.

`nbfb` treats a t-BFB cycle as a fold-back t-BFB cycle when the last segment of
`s` and the first segment of `x` lie on one chromosome in opposite orientations,
with breakpoints at most 50 kb apart (`--foldback-span`).

### NBFB Walk

A **breakpoint graph** has a segment edge per segment, reference edges between
segments that are adjacent in the reference genome, and junctions. The
junctions that join two left ends or two right ends close together on one
chromosome are the **fold-back junctions**. Every edge has a copy number.

An **alternating walk** alternates segment edges with reference edges or
junctions and traverses every edge at most its copy number of times. The walk
spells its target. A walk from a node `v` is **maximal** if it maximizes the
number of traversals of fold-back junctions, and then the number of traversals
of segment edges, among the alternating walks from `v`.

Given a breakpoint graph and a node `v`, **NBFB Walk** seeks a maximal walk from
`v` and a nested history for its target that minimizes the cost of
Interchromosomal NBFB Distance.

## Workflow

### From an assembled target

Write the sources and the target as an input file (see
[Sequence format](#sequence-format)) and run:

```bash
cpp/nbfb -i instance.txt
```

The last line of the output is the minimum cost of a history, an integer
(`total`) or `inf`. The lines before it give `ft` (`foldback_t`), `t`
(`t_cycles`), the cycles of an optimal history and the number of optimal
histories.

### From a breakpoint graph

1. Build the breakpoint graph of a locus with `tracebfb graph`, or write one in
   the format of [`docs/file_formats.md`](docs/file_formats.md).
2. Count the maximal walks from the start node `v`, for example the left end
   `1.L` of segment 1:

   ```bash
   cpp/lwalk graph.tsv --starts 1.L --max-foldbacks
   ```

   List them all, or sample them uniformly when there are too many:

   ```bash
   cpp/lwalk graph.tsv --starts 1.L --max-foldbacks --enumerate -o walks
   cpp/lwalk graph.tsv --starts 1.L --max-foldbacks --sample 10000 -o walks
   ```

   `walks.walks.txt` gets one walk per line as signed segment ids, and
   `walks.targets.txt` the target of every walk.
3. Solve Interchromosomal NBFB Distance for the target of every walk, with the
   strings of the windows as source lines in `sources.txt`:

   ```bash
   cpp/nbfb --sources sources.txt --targets walks.targets.txt -o nbfb.tsv
   ```

   `nbfb.tsv` gets one row per walk with its cost (`foldback_t`, `t_cycles`,
   `total`) and an optimal history. The walks of minimum cost solve NBFB Walk.

See [`cpp/README.md`](cpp/README.md) for every option of `nbfb` and `lwalk`,
the output formats and the algorithms.

## Building a breakpoint graph

```bash
tracebfb graph \
    --region chr4:127186850-129234631 \
    --junctions tumor_junctions.tsv \
    --normal-junctions normal_junctions.tsv \
    --bam tumor.bam \
    --normal-bam normal.bam \
    --haploid-depth 15 \
    --cn wakhan_hp1.bed --cn wakhan_hp2.bed \
    --cytoband cytoBand.hg38.txt.gz \
    --name HCC1937_chr4 \
    --out graph_dir
```

`tracebfb graph` builds the breakpoint graph of one locus and writes
`graph_dir/graph.tsv`:

- **Windows.** `--region` gives one window of the locus and may be repeated.
- **Junctions.** The junction files hold the junctions called from split reads,
  as TSV with a header and the columns `chrA`, `posA`, `sideA`, `chrB`, `posB`,
  `sideB` and `support` (the read support). The side is `R` when the segment
  lies to the left of the position. A tumor junction is kept if its read
  support is at least 5, it lies more than 1 Mb from a centromere (with
  `--cytoband`), it is not a deletion of at most 600 bp, and no normal junction
  of the same form lies within 200 bp. A junction is a fold-back junction when
  its positions are at most 50 kb apart, unless an opposite junction within
  10 kb indicates an inversion.
- **Segments.** The windows are cut at the positions of the junctions.
- **Copy numbers.** `--haploid-depth` is the haploid read depth, the tumor read
  depth of one copy. `--cn` gives copy-number calls (BED, copy number in column
  5), such as the two haplotypes of Wakhan, summed over files. The copy numbers
  minimize a depth error, a junction deviance and a prior on the read depth per
  copy, make the graph semi-balanced, and give every junction a copy number of
  at least 1.

`--reads SET=BAM` adds a read set whose subpath constraints, the segments that
every read crosses, are written to `graph_dir/subpaths.tsv`. With
`--use-subpaths`, the subpath constraints are lower bounds in the fit of the
copy numbers.

The graph and subpath formats are described in
[`docs/file_formats.md`](docs/file_formats.md).

## Sequence format

A sequence is a semicolon-separated list of genomic segments:

```text
chr1:100-200[R];chr1:200-300[R]
```

Each segment has the form `<chromosome>:<start>-<end>[R|C]`, where `R` is
the reference orientation and `C` the reverse-complement orientation. The
reverse complement of `chr1:100-200[R]` is `chr1:200-100[C]`.

An input file lists the sources followed by the target, one non-empty sequence
per line:

```text
chr1:100-200[R];chr1:200-300[R]
chr2:500-600[R]
chr1:100-200[R];chr1:200-300[R];chr1:300-200[C];chr2:500-600[R]
```

`nbfb -i` and `tracebfb distance` read the same input files.

## Histories that need not be nested

```bash
tracebfb distance --input sequences.txt --out result --colored
```

`tracebfb distance` computes the minimum cost of a history over the sources as
listed, nested or not, by backward search from the target (exponential time),
and prints the total number of cycles. The cost is `(t, total)`, as with
`nbfb --general --no-foldback-first`. A nested history can cost more than a
history that is not nested, but only with overhangs of unrealistic length.

The command creates:

- `result.graph`: the BFB graph of all optimal histories, in the text format
  below.
- `result.dot`: the BFB graph in Graphviz DOT format, edges labeled by cycle.
- `result.paths`: every optimal history.
- `result.paths.ansi`: the same histories with ANSI coloring when
  `--colored` is supplied.

Within a path file, a history of `k` cycles occupies `k` lines: the
string produced by each of its first `k - 1` cycles, then the target. The
empty start string is omitted. A blank line separates alternative
histories. A target without a history produces an empty path file.

ANSI path files color each string by the cycle applied to it: `p` red,
`q` green, and `r` purple for p-, i-, and f-BFB cycles, and the whole
string cyan for a t-BFB cycle. Removing the ANSI escape sequences produces
the corresponding uncolored path file.

### BFB graph text format

The first line is the target and the second line the tab-separated sources:

```text
<target>
<source_0>\t<source_1>...
```

Each remaining line is an adjacency row:

```text
<node>\t<child>,<kind>,<x>,<y>\t<child>,t,<source>...
```

A node is written as the set of appended source indices in braces followed
by the string, for example `{0,2}chr1:100-200[R]`; the start state is
`{}`. `x` and `y` encode the decomposition `p = node[:x]`, `q = node[x:y]`,
and `r = node[y:]`.

The end state has no adjacency row. A graph of a target that is empty or has
no history contains only the `{}` row without an edge.

## Simulating histories

```bash
tracebfb simulate \
    --seed 11 \
    --n-sources 2 \
    --n-segments 5 \
    --cycles 6 \
    --short-deletions \
    --out simulation \
    --colored
```

The sources are generated on distinct chromosomes. The history opens with a
t-BFB cycle, appends every source once at random positions, and fills the
remaining cycles with random p-, i- and f-BFB cycles. With
`--short-deletions`, every overhang of an f-BFB cycle is shorter than the
string that the preceding cycle added, so the history is nested.

The command creates:

- `simulation.input`: the sources, then the target, ready for
  `tracebfb distance` and `nbfb -i`.
- `simulation.paths`: the string after every cycle, one per line.
- `simulation.paths.ansi`: the same strings with ANSI coloring when
  `--colored` is supplied.

Other options are `--min-len` and `--max-len` for segment lengths.

## Python API

```python
from tracebfb import BFB_dist
from tracebfb.sequence import Sequence

sources = [
    Sequence.from_string("chr1:100-200[R];chr1:200-300[R]"),
    Sequence.from_string("chr2:500-600[R]"),
]
target = sources[0].bfb_target(0, 1, False) + sources[1]

graph = BFB_dist(sources, target)
print(graph.calculate_bfb_distance())  # 3: t, i, t

for path in graph.get_paths():
    for state, child in zip(path, path[1:]):
        operation = graph.adjacency.get(state, {}).get(child)
        print(state.sequence, state.translocations, operation)

with open("result.graph", "w", encoding="utf-8") as output:
    graph.write_text(output)
```

`BFB_dist` takes a list of sources and a target and returns a `BFBGraph`.
Graph nodes are `State` objects pairing a string with the number of t-BFB
cycles applied so far (`translocations`). `graph.start` is the empty state, and
`graph.ends` is the set of states that spell the target. Edges carry an
`Operation` with a `kind` in `p`, `i`, `f`, `t`, the decomposition indices
`x` and `y` for p-, i- and f-BFB cycles, and the `source` index for t-BFB
cycles.

Read and validate a saved BFB graph with:

```python
from tracebfb import BFBGraph

with open("result.graph", encoding="utf-8") as input_file:
    graph = BFBGraph.from_text(input_file)
```

`BFBGraph.from_text` validates the graph because the file is treated as
untrusted input. It verifies the endpoints, adjacency structure, operations,
reachability, terminal behavior, and equal start-to-end path lengths.

`BFBGraph` caches its distance and paths. Treat its sources, target,
adjacency, and returned paths as immutable.

Other modules:

- `tracebfb.foldback_first.foldback_first_BFB_dist(sources, target)`: the
  minimum cost `(ft, t, total)` of a history that need not be nested, with one
  optimal history. `cpp/nbfb --general` is the nested counterpart.
- `tracebfb.breakpoint_graph.BreakpointGraph`: the breakpoint graph of a
  target over the reference edges of the sources (`from_target`), written in
  the graph format (`to_file`) or the input format of `cpp/lwalk`
  (`write_lwalk`).
- `tracebfb.graph_file.GraphFile`: reads, writes and validates the graph
  format.
- `tracebfb.graph_construction.build_graph`: the function behind
  `tracebfb graph`.

## Tests

Run the test suite from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src \
python -m unittest discover -s tests -v
```

The tests check the backward searches on simulated histories, compare `nbfb`
with a brute-force search that applies the definition of a nested history
literally, and check `lwalk` on breakpoint graphs of known targets. The tests
of `nbfb` and `lwalk` build them with `make` and are skipped when the build
fails.
