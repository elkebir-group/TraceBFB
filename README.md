# BFB-dist

BFB-dist simulates breakage-fusion-bridge (BFB) histories over a set of
source sequences and infers minimum-length histories that produce a target
sequence.

## Installation

BFB-dist requires Python 3.7 or newer and has no runtime dependencies.

```bash
python -m pip install .
```

This installs the `bfb-dist` command and the `bfb_dist` Python package.

## Model

A history starts from the empty sequence and applies BFB cycles. A source
may be appended any number of times. For a
sequence `P + Q + R` with nonempty `Q`, the cycles are:

| kind | name                          | result                                        | constraint    |
|------|-------------------------------|-----------------------------------------------|---------------|
| `p`  | proper BFB                    | `P + Q + reverse_complement(Q)`               | `R` empty     |
| `i`  | inverted-part deletion        | `P + Q + R + reverse_complement(Q)`           | `R` nonempty  |
| `f`  | forward-part deletion         | `P + Q + reverse_complement(R) + reverse_complement(Q)` | `R` nonempty |
| `t`  | translocation                 | `sequence + source`                           |               |

The **BFB distance** of a set of sources and a target is the cost of an
optimal history producing the target: minimum first in the number of `t`
cycles and then in the total number of cycles. It is `0` only for an empty
target, and infinite when no history exists.

A history is **nested** if every sequence it produces is `t[1, psi] y`, a
prefix of the target followed by a string `y` that is empty, the reverse
complement of a substring of the target, or a source. The **NBFB distance**
is the cost of an optimal nested history.

## Sequence format

A sequence is a semicolon-separated list of genomic segments:

```text
chr1:100-200[R];chr1:200-300[R]
```

Each segment has the form `<chromosome>:<start>-<end>[R|C]`, where `R` is
the reference orientation and `C` the reverse-complement orientation.

## Command-line use

### Simulate a history

```bash
bfb-dist simulate \
    --seed 11 \
    --n-sources 2 \
    --n-segments 5 \
    --cycles 6 \
    --short-deletions \
    --out simulation \
    --colored
```

The sources are generated on distinct chromosomes. The history opens with a
`t` cycle, appends every source once at random positions, and fills the
remaining cycles with random `p`, `i`, and `f` cycles. With
`--short-deletions`, every `f` cycle removes fewer segments than the
preceding cycle added, which guarantees a nested history.

The command creates:

- `simulation.input`: the sources followed by the target, one per line,
  ready for `bfb-dist distance`.
- `simulation.paths`: the sequence after every cycle, one per line.
- `simulation.paths.ansi`: the same sequences with ANSI coloring when
  `--colored` is supplied.

Other options are `--min-len` and `--max-len` for segment lengths.

### Calculate BFB distance

Create an input file whose lines are the sources followed by the target,
one nonempty sequence per line:

```text
chr1:100-200[R];chr1:200-300[R]
chr2:500-600[R]
chr1:100-200[R];chr1:200-300[R];chr1:300-200[C];chr2:500-600[R]
```

Then run:

```bash
bfb-dist distance --input sequences.txt --out result --colored
```

The distance is printed to standard output. Pass `--mode` to select the
history class and algorithm:

- `--mode bfb` (default): the general BFB distance, by backward search from
  the target. Exponential time.
- `--mode nbfb`: the NBFB distance, by breadth-first search over the NBFB
  graph. Returns every optimal nested history. For large instances use the
  C++ tool in `cpp/` (`nbfb`), which runs in `O(n^4)` time when the sources
  are all blocks of a few strings.

The command creates:

- `result.graph`: the BFB graph in its text adjacency-list format.
- `result.dot`: the graph in Graphviz DOT format, edges labeled by kind.
- `result.paths`: every minimum-length history in the graph.
- `result.paths.ansi`: the same histories with ANSI coloring when
  `--colored` is supplied.

Within a path file, a history of length `k` occupies `k` lines: the
sequence produced by each of its first `k - 1` cycles, then the target. The
empty start sequence is omitted. A blank line separates alternative
histories. An unreachable target produces an empty path file.

ANSI path files color each sequence by the cycle applied to it: `P` red,
`Q` green, and `R` purple for `p`, `i`, and `f` cycles, and the whole
sequence cyan for a `t` cycle. Removing the ANSI escape sequences produces
the corresponding uncolored path file.

## Python API

```python
from bfb_dist import BFB_dist, NBFB_dist
from bfb_dist.sequence import Sequence

sources = [
    Sequence.from_string("chr1:100-200[R];chr1:200-300[R]"),
    Sequence.from_string("chr2:500-600[R]"),
]
target = sources[0].bfb_target(0, 1, False) + sources[1]

graph = BFB_dist(sources, target)
print(graph.calculate_bfb_distance())  # 3: t, p, t

for path in NBFB_dist(sources, target).get_paths():
    for state, child in zip(path, path[1:]):
        operation = graph.adjacency.get(state, {}).get(child)
        print(state.sequence, state.translocations, operation)

with open("result.graph", "w", encoding="utf-8") as output:
    graph.write_text(output)
```

Both functions take a list of sources and a target and return a
`BFBGraph`. Graph nodes are `State` objects pairing a sequence with the
number of t-BFB cycles applied so far (`translocations`). `graph.start` is
the empty state, and `graph.ends()` lists the states that spell the target.
Edges carry an
`Operation` with a `kind` in `p`, `i`, `f`, `t`, the decomposition indices
`x` and `y` for fold-back kinds, and the `source` index for `t`.

Read and validate a saved graph with:

```python
from bfb_dist import BFBGraph

with open("result.graph", encoding="utf-8") as input_file:
    graph = BFBGraph.from_text(input_file)
```

`BFBGraph.from_text` validates the graph because the file is treated as
untrusted input. It verifies the endpoints, adjacency structure, operations,
reachability, terminal behavior, and equal start-to-end path lengths.

`BFBGraph` caches its distance and paths. Treat its sources, target,
adjacency, and returned paths as immutable.

## Graph text format

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
by the sequence, for example `{0,2}chr1:100-200[R]`; the start state is
`{}`. `x` and `y` encode the decomposition `P = node[:x]`, `Q = node[x:y]`,
and `R = node[y:]`.

The end state has no adjacency row. A zero-distance or unreachable graph
contains only the `{}` row without an edge.

## Tests

Run the test suite from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src \
python -m unittest discover -s tests -v
```

The suite checks the nested and general searches against each other on
simulated histories and checks that short-deletion histories are nested.
