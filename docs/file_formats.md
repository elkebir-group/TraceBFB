# File formats — DRAFT

`bfbdist` reads and writes three files:

| File | Extension | Holds | Needed by |
|---|---|---|---|
| Breakpoint graph | `.tsv` | segments, reference edges, junctions, free ends | walk enumeration (`lwalk`), graph distances |
| Sources | `.txt` | the source strings `g_1, ..., g_r` | NBFB and BFB distances |
| Subpath constraints | `.tsv` | the oriented segments that reads cross, with read counts | optional analyses only |

The graph never depends on the other two files. All three identify segments the same way.

## Segment ids and strings

A **segment id** is a positive integer (`1`, `2`, …). The negated id denotes the reverse
complement (`-2`).

A **string** is a sequence of segment ids and reverse complements, separated by single spaces:
`1 2 3 -3 -2`. The same notation is used for sources, targets, walks and subpath constraints.

Segment `X` has two **ends**, `X.L` and `X.R` (the paper's left node `v⁻` and right node `v⁺`).
Traversing `X` from `X.L` to `X.R` spells `X`; traversing it from `X.R` to `X.L` spells `-X`.

---

## 1. Breakpoint graph (`.tsv`)

The breakpoint graph of the paper (`def:bpgraph`): segments with copy numbers, reference edges,
junctions (variant edges), and free ends. Nothing else.

### Records

Tab-separated, one record per line, no column header. The first field is the record type. Lines
starting with `#` are comments. Optional fields come last as `key=value`.

| Type | Fields | Meaning |
|---|---|---|
| `H` | `key=value ...` | Header, first record, `version=1`. Any other `key=value` (`sample`, `locus`, `reference`, ...) is kept as metadata. |
| `S` | `id chrom start end strand copies` | Segment: the sequence `chrom:start-end` (0-based half-open), read on `strand` (`+` or `-`), with `copies` ≥ 0. `X.L` is where that reading starts, `X.R` where it ends. |
| `R` | `end1 end2 multiplicity` | Reference edge. |
| `J` | `name end1 end2 multiplicity foldback` | Junction. `foldback` is `0` or `1`. |
| `F` | `end count` | Free end: `count` copies of segment end `end` are joined to no edge (molecule ends). Where a walk starts is a parameter of the walk search, not part of the graph. |

### Validity

A reader rejects the file unless:

1. Segment ids are unique, junction names are unique, and every `R`, `J`, `F` refers to
   declared segments.
2. Every segment end has at most one `R`.
3. Multiplicities: `R` ≥ 1, `J` ≥ 1, `F` count ≥ 1.
4. **Balance** at every end `v` of every segment `X`:
   `copies(X)` = Σ multiplicities of the edges at `v` (a junction `v`–`v` counts twice) + free ends at `v`.
5. The total number of free ends is even and ≥ 2.

### Derived, not stored

- **Class.** `single` if exactly 2 free ends (one molecule), `multi` otherwise.
- **Walks.** A walk starts at a free end, alternates segment traversals and edges, and uses each
  segment at most `copies` times and each edge at most its multiplicity. It spells a string.

### Example: HCC1937 chr4q

```
H	version=1	sample=HCC1937	locus=HCC1937_chr4_128.19	reference=GRCh38
S	1	chr4	127186850	128186850	+	2
S	2	chr4	128186850	128187424	+	3
S	3	chr4	128187424	128232719	+	4
S	4	chr4	128232719	128234631	+	2
S	5	chr4	128234631	129234631	+	0
R	1.R	2.L	2
R	2.R	3.L	3
R	3.R	4.L	2
J	j0	2.L	3.L	1	1
J	j1	3.R	4.R	2	1
F	1.L	2
```

Balance at `3.R`: `copies(3)` = 4 = `R 3.R 4.L` (2) + `J j1` (2).

---

## 2. Sources (`.txt`)

The source strings `g_1, ..., g_r` of a distance computation (`\mathcal{S}` in the paper).

- One source per line, as a string (signed segment ids, negative = reverse complement). Lines starting with `#` are comments. Blank lines are ignored.
- A source may contain each segment at most once.
- When used with a graph, every id must be a segment of the graph.

Example (HCC1937 chr4q, one source):

```
1 2 3 4 5
```

Example (two sources, the second read on the reverse strand):

```
6 7 8
-11 -10 -9
```

---

## 3. Subpath constraints (`.tsv`)

The oriented segments that reads cross, one row per distinct string. Optional input: only
analyses that score walks against reads use it.

Tab-separated, with a column header. Lines starting with `#` are comments.

| Column | Required | Meaning |
|---|---|---|
| `string` | yes | the segments one read crosses, in order (a string) |
| `reads` | yes | number of reads that cross exactly this string |
| `read_set` | no | the read set (`HiFi`, `ONT_R10`, `ONT_E821`, ...) |
| `example_read` | no | the name of one such read |

- A string and its reverse complement are the same constraint. Each is written once, in the
  reading with more forward segments (ties: the smaller reading as integers), e.g. `1 2 3`, not
  `-3 -2 -1`.
- With `read_set`, a string may appear once per read set.
- When used with a graph, every id must be a segment of the graph. Whether consecutive segments
  are joined by an edge is checked by the analysis, not by the reader.

Example (HCC1937 chr4q):

```
string	reads	read_set
3 4 -3 2 3 -4 -3	2	ONT_E821
```

---

## Implementation in `bfbdist`

- `bfbdist.graph_file.GraphFile`: `read(path)`, `write(path)`, `validate()` for the graph format.
- `bfbdist graph` (`bfbdist.graph_construction`) builds a graph from long reads and writes
  `graph.tsv` and `subpaths.tsv`.
- `bfbdist.breakpoint_graph.BreakpointGraph.to_file()` gives the graph of an assembled target.
- `cpp/lwalk` reads the graph format (and still the older `segments/bonds/starts` format).
- Not yet: `cpp/nbfb` reading sources in this format.
