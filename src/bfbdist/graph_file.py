"""The breakpoint graph file (`docs/file_formats.md`): read, write and validate.

A graph holds segments with copy numbers, reference edges, junctions and free
ends. Segments have positive integer ids. A segment end is written ``X.L`` (where
the reading of segment ``X`` starts) or ``X.R`` (where it ends).
"""

from collections import Counter
from dataclasses import dataclass, field

VERSION = "1"

# The fold-back geometry of the junction screen, and the default of `nbfb --foldback-span`.
FOLDBACK_SPAN = 50_000


def segment_ids(n):
    """1, 2, ..., n: the segment ids. The reverse complement of segment k is -k."""
    return [str(k) for k in range(1, n + 1)]


def is_foldback(a, b, span=FOLDBACK_SPAN):
    """Whether a junction between ends a and b, each (chrom, position, kind), is a fold-back:
    two ends of the same kind on one chromosome at most span bp apart."""
    return a[2] == b[2] and a[0] == b[0] and abs(a[1] - b[1]) <= span


@dataclass
class GraphSegment:
    name: str
    chrom: str
    start: int
    end: int
    strand: str
    copies: int

    def token(self):
        """The segment as a `chrom:start-end[R]` token of the C++ tools."""
        if self.chrom == "*":
            return self.name
        if self.strand == "+":
            return f"{self.chrom}:{self.start}-{self.end}[R]"
        return f"{self.chrom}:{self.end}-{self.start}[C]"


@dataclass
class GraphFile:
    header: dict = field(default_factory=dict)
    segments: list = field(default_factory=list)   # GraphSegment
    reference: list = field(default_factory=list)  # (end1, end2, multiplicity)
    junctions: list = field(default_factory=list)  # (name, end1, end2, multiplicity, foldback)
    free_ends: list = field(default_factory=list)  # (end, count)

    def write(self, path):
        lines = ["\t".join(["H", f"version={VERSION}"]
                           + [f"{k}={v}" for k, v in self.header.items() if k != "version"])]
        for s in self.segments:
            coords = ["*", "*", "*"] if s.chrom == "*" else [s.chrom, str(s.start), str(s.end)]
            lines.append("\t".join(["S", s.name, *coords, s.strand, str(s.copies)]))
        for u, v, m in self.reference:
            lines.append("\t".join(["R", u, v, str(m)]))
        for name, u, v, m, fb in self.junctions:
            lines.append("\t".join(["J", name, u, v, str(m), str(int(fb))]))
        for end, count in self.free_ends:
            lines.append("\t".join(["F", end, str(count)]))
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")

    @classmethod
    def read(cls, path):
        graph = cls()
        with open(path, encoding="utf-8") as handle:
            for number, raw in enumerate(handle, 1):
                line = raw.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                kind, *f = line.split("\t")
                try:
                    if kind == "H":
                        graph.header = dict(x.split("=", 1) for x in f)
                    elif kind == "S":
                        name, chrom, start, end, strand, copies = f[:6]
                        if chrom == "*":
                            start = end = 0
                        graph.segments.append(GraphSegment(
                            name, chrom, int(start), int(end), strand, int(copies)))
                    elif kind == "R":
                        graph.reference.append((f[0], f[1], int(f[2])))
                    elif kind == "J":
                        graph.junctions.append((f[0], f[1], f[2], int(f[3]), int(f[4])))
                    elif kind == "F":
                        graph.free_ends.append((f[0], int(f[1])))
                    else:
                        raise ValueError(f"unknown record type {kind!r}")
                except (ValueError, IndexError) as error:
                    raise ValueError(f"{path}:{number}: {error}") from None
        graph.validate()
        return graph

    def validate(self):
        """Raise ValueError unless the graph satisfies the rules of the format."""
        if self.header.get("version", VERSION) != VERSION:
            raise ValueError("header must give version=1")
        copies = {}
        for s in self.segments:
            if s.name in copies:
                raise ValueError(f"segment {s.name} declared twice")
            if not s.name.isdigit() or s.name.startswith("0"):
                raise ValueError(f"segment id {s.name} must be a positive integer")
            if s.strand not in "+-" or len(s.strand) != 1:
                raise ValueError(f"segment {s.name}: strand must be + or -")
            if s.copies < 0:
                raise ValueError(f"segment {s.name}: negative copies")
            copies[s.name] = s.copies

        def check(end):
            name, _, side = end.rpartition(".")
            if name not in copies or side not in ("L", "R"):
                raise ValueError(f"unknown segment end {end}")

        degree = Counter()
        seen = set()
        for u, v, m in self.reference:
            check(u)
            check(v)
            for x in (u, v):
                if x in seen:
                    raise ValueError(f"segment end {x} has two reference edges")
                seen.add(x)
            if m < 1:
                raise ValueError(f"reference edge {u} {v}: multiplicity must be >= 1")
            degree[u] += m
            degree[v] += m
        names = set()
        for name, u, v, m, fb in self.junctions:
            check(u)
            check(v)
            if name in names:
                raise ValueError(f"junction {name} declared twice")
            names.add(name)
            if m < 1:
                raise ValueError(f"junction {name}: multiplicity must be >= 1")
            if fb not in (0, 1):
                raise ValueError(f"junction {name}: foldback must be 0 or 1")
            degree[u] += m
            degree[v] += m
        total = 0
        for end, count in self.free_ends:
            check(end)
            if count < 1:
                raise ValueError(f"free end {end}: count must be >= 1")
            degree[end] += count
            total += count
        for name, c in copies.items():
            for side in ("L", "R"):
                if degree[f"{name}.{side}"] != c:
                    raise ValueError(
                        f"unbalanced at {name}.{side}: {c} copies, "
                        f"{degree[f'{name}.{side}']} edge and free-end copies")
        if total < 2 or total % 2:
            raise ValueError(f"{total} free ends: must be even and at least 2")
