# BFB-dist

BFB-dist simulates breakage-fusion-bridge (BFB) cycles and infers every
minimum-distance BFB path between two sequences.

## Installation

BFB-dist requires Python 3.7 or newer and has no runtime dependencies.

```bash
python -m pip install .
```

This installs the `bfb-dist` command and the `bfb_dist` Python package.

## Sequence format

A sequence is a semicolon-separated list of genomic segments:

```text
chr1:100-200[R];chr1:200-300[R]
```

Each segment has the form:

```text
<chromosome>:<start>-<end>[R|C]
```

`R` is the reference orientation and `C` is the reverse-complement
orientation.

For a source sequence `P + Q + R`, one BFB operation produces either:

```text
P + Q + R + reverse_complement(Q)
```

or, for a forward deletion:

```text
P + Q + reverse_complement(R) + reverse_complement(Q)
```

## Command-line use

### Simulate BFB cycles

```bash
bfb-dist simulate \
    --seed 11 \
    --n-segments 5 \
    --cycles 3 \
    --out simulation \
    --colored
```

The command creates:

- `simulation.paths`: the reference and every successive simulated sequence,
  one complete sequence per line.
- `simulation.paths.ansi`: the same sequences with ANSI coloring when
  `--colored` is supplied.

Other simulation options are:

```text
--min-len INTEGER
--max-len INTEGER
```

### Calculate BFB distance

Create an input file containing the source and target on two consecutive
nonempty lines:

```text
chr1:100-200[R];chr1:200-300[R]
chr1:100-200[R];chr1:200-300[R];chr1:300-200[C]
```

Then run:

```bash
bfb-dist distance --input sequences.txt --out result --colored
```

By default this computes the general, unrestricted BFB distance. Pass
`--mode` to restrict the search to a more biologically-focused variant:

- `--mode bfb` (default): every p/i/f-BFB history reaching the target.
- `--mode nbfb`: only *nested* histories, where every forward-part deletion
  stays within the material the immediately preceding cycle created.
- `--mode nbfb-star`: only *strongly nested* histories, a polynomial-time
  restriction of `nbfb` to histories whose intermediate sequences are always
  foldable back onto the target. `nbfb-star` distance is always at least the
  `nbfb` distance, and can be unreachable (infinite) even when `nbfb` and
  `bfb` are both finite.

The command creates:

- `result.graph`: the BFB graph in its text adjacency-list format.
- `result.dot`: the graph in Graphviz DOT format.
- `result.paths`: every minimum-distance path.
- `result.paths.ansi`: the same paths with ANSI coloring when `--colored` is
  supplied.

Within a path file, every nonempty line is one complete sequence. A blank line
separates alternative source-to-target paths. A zero-distance graph contains
one line with the source. An unreachable graph produces an empty path file.

ANSI path files color the P, Q, and R portions of each source or intermediate
sequence. Removing the ANSI escape sequences produces the corresponding
uncolored path file.

## Python API

`BFB_dist`, `NBFB_dist`, `NBFB_star_dist`, and `BFBGraph` are available
directly from the package. `NBFB_dist` and `NBFB_star_dist` share `BFB_dist`'s
signature and return the same `BFBGraph` type, restricted to nested and
strongly nested histories respectively:

```python
from bfb_dist import BFB_dist, BFBGraph, NBFB_dist, NBFB_star_dist
from bfb_dist.sequence import Sequence

source = Sequence.from_string(
    "chr1:100-200[R];chr1:200-300[R]"
)
target = source.bfb_target(0, 1, False)

graph = BFB_dist(source, target)
print(graph.calculate_bfb_distance())

nested_graph = NBFB_dist(source, target)
print(nested_graph.calculate_bfb_distance())

strongly_nested_graph = NBFB_star_dist(source, target)
print(strongly_nested_graph.calculate_bfb_distance())

with open("result.graph", "w", encoding="utf-8") as output:
    graph.write_text(output)

with open("result.dot", "w", encoding="utf-8") as output:
    graph.write_dot(output)

with open("result.paths", "w", encoding="utf-8") as output:
    output.write(graph.write_paths(colored=False))
```

Read and validate a saved graph with:

```python
from bfb_dist import BFBGraph

with open("result.graph", encoding="utf-8") as input_file:
    graph = BFBGraph.from_text(input_file)
```

`BFBGraph.from_text` validates the graph because the file is treated as
untrusted input. It verifies the endpoints, adjacency structure, BFB
operations, reachability, terminal behavior, and equal source-to-target path
lengths.

Graphs returned by `BFB_dist` are constructed from minimum-distance layers, so
they satisfy these invariants without running the file-loader validation.

`BFBGraph` caches its distance and paths. Treat its source, target, adjacency,
and returned paths as immutable.

## Graph text format

The first line contains the source and target separated by a tab:

```text
<source>\t<target>
```

Each remaining line is an adjacency row:

```text
<node>\t<child>,<x>,<y>,<is_forward_deletion>...
```

`x` and `y` encode the decomposition `P = node[:x]`, `Q = node[x:y]`, and
`R = node[y:]`. `is_forward_deletion` is `0` or `1`.

The target has no adjacency row. A zero-distance or unreachable graph contains
only the source row without an edge.

## Tests

Run the test suite from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src \
python -m unittest discover -s tests -v
```
