"""Breakpoint graphs and their walks.

A node is a segment end: ``(segment, "h")`` is the start and ``(segment, "t")``
the end of the segment in its forward orientation. A reference edge joins two
ends that are adjacent in a source. A variant edge (a junction) joins two ends
that are adjacent in the target but in no source; a fold-back whose two ends
coincide is a self-loop. Every edge has a multiplicity and every segment a
copy number.

A walk enters a segment at one end, traverses the segment, and leaves at the
other end, so the walk spells a sequence. ``cpp/lwalk`` enumerates, counts,
and samples the longest walks of a graph written by
:meth:`BreakpointGraph.write_lwalk`.
"""

from collections import Counter

from bfbdist.graph_file import FOLDBACK_SPAN, GraphFile, GraphSegment, is_foldback, segment_ids


def _is_forward(segment):
    if hasattr(segment, "complement"):
        return not segment.complement
    return segment.letter.isupper()


def _forward(segment):
    return segment if _is_forward(segment) else segment.reverse_complement()


def enter(segment):
    """Return the node through which a walk enters ``segment``."""
    if _is_forward(segment):
        return (segment, "h")
    return (segment.reverse_complement(), "t")


def leave(segment):
    """Return the node through which a walk leaves ``segment``."""
    if _is_forward(segment):
        return (segment, "t")
    return (segment.reverse_complement(), "h")


def _edge(u, v):
    return frozenset((u, v))


def _other(edge, node):
    if len(edge) == 1:
        return node
    return next(x for x in edge if x != node)


def _adjacencies(sequence):
    return [
        _edge(leave(a), enter(b))
        for a, b in zip(sequence.segments, sequence.segments[1:])
    ]


class BreakpointGraph:
    """A breakpoint graph with two free ends.

    Attributes:
        copy_numbers: ``Counter`` of forward segments.
        junctions: ``Counter`` of variant edges, each a frozenset of nodes.
        reference: ``Counter`` of reference edges.
        start: the node at which every walk begins.
        end: the node at which every walk ends.
    """

    def __init__(self, copy_numbers, junctions, reference, start, end):
        self.copy_numbers = Counter(copy_numbers)
        self.junctions = Counter(junctions)
        self.reference = Counter(reference)
        self.start = start
        self.end = end

    @staticmethod
    def degree(node, edges):
        """Return the edge ends at ``node``, with a self-loop counted twice."""
        return sum(
            m * (2 if len(e) == 1 else 1) for e, m in edges.items() if node in e
        )

    @classmethod
    def from_target(cls, sources, target):
        """Return the graph that junction-spanning reads and read depth report
        for ``target``.

        The copy numbers, the junctions with their multiplicities, and the
        two free ends are read from the target. The multiplicity of each
        reference edge follows from balance at its two ends.

        Raises:
            ValueError: If the balance at the two ends of a reference edge
                disagrees, or disagrees with the target.
        """
        reference_edges = {e for s in sources for e in _adjacencies(s)}
        copy_numbers = Counter(_forward(s) for s in target.segments)
        adjacency = Counter(_adjacencies(target))
        junctions = Counter(
            {e: m for e, m in adjacency.items() if e not in reference_edges}
        )
        start, end = enter(target[0]), leave(target[-1])
        free = Counter([start, end])
        reference = Counter()
        for e in reference_edges:
            if len(e) == 1:
                continue
            u, v = tuple(e)
            mu = copy_numbers[u[0]] - free[u] - cls.degree(u, junctions)
            mv = copy_numbers[v[0]] - free[v] - cls.degree(v, junctions)
            if mu != mv or mu < 0:
                raise ValueError(f"unbalanced reference edge {set(e)}: {mu} vs {mv}")
            if mu != adjacency.get(e, 0):
                raise ValueError(
                    f"balance gives {mu} for {set(e)}, the target "
                    f"{adjacency.get(e, 0)}"
                )
            if mu:
                reference[e] = mu
        return cls(copy_numbers, junctions, reference, start, end)

    def to_file(self, names=None, foldback_span=FOLDBACK_SPAN, header=None):
        """Return the graph as a :class:`bfbdist.graph_file.GraphFile`.

        ``names`` maps a forward segment to its id, by default ``1``, ``2``,
        ... in order of first occurrence. A segment without coordinates (a
        letter) is written with ``*`` coordinates.
        """
        segments = list(self.copy_numbers)
        if names is None:
            names = dict(zip(segments, segment_ids(len(segments))))

        def end(v):
            return f"{names[v[0]]}.{'L' if v[1] == 'h' else 'R'}"

        graph = GraphFile(header=dict(header or {}))
        for s in segments:
            if hasattr(s, "chr"):
                graph.segments.append(GraphSegment(
                    names[s], s.chr, min(s.start, s.end), max(s.start, s.end),
                    "-" if s.start > s.end else "+", self.copy_numbers[s]))
            else:
                graph.segments.append(GraphSegment(names[s], "*", 0, 0, "+", self.copy_numbers[s]))
        for e, m in self.reference.items():
            u = next(iter(e))
            graph.reference.append((end(u), end(_other(e, u)), m))
        for k, (e, m) in enumerate(self.junctions.items()):
            u = next(iter(e))
            graph.junctions.append((f"j{k}", end(u), end(_other(e, u)), m,
                                    int(_foldback(e, foldback_span))))
        for v, count in Counter([self.start, self.end]).items():
            graph.free_ends.append((end(v), count))
        graph.validate()
        return graph

    def write_lwalk(self, path, names=None, foldback_span=FOLDBACK_SPAN):
        """Write the graph in the input format of ``cpp/lwalk``.

        ``names`` maps a forward segment to its name in the walk letters, by
        default ``s0``, ``s1``, ... in order of first occurrence. A node is
        ``2 * i`` for the start and ``2 * i + 1`` for the end of segment ``i``.
        A junction is flagged as a fold-back when it joins two ends of the
        same kind on one chromosome at most ``foldback_span`` bp apart.
        """
        segments = list(self.copy_numbers)
        index = {s: i for i, s in enumerate(segments)}
        if names is None:
            names = {s: f"s{i}" for i, s in enumerate(segments)}

        def node(v):
            return 2 * index[v[0]] + (0 if v[1] == "h" else 1)

        lines = [f"segments {len(segments)}"]
        lines += [f"{names[s]} {self.copy_numbers[s]} {s}" for s in segments]
        bonds = []
        for edges, kind in ((self.reference, "ref"), (self.junctions, "junction")):
            for k, (e, m) in enumerate(edges.items()):
                u = next(iter(e))
                v = _other(e, u)
                fb = int(kind == "junction" and _foldback(e, foldback_span))
                bonds.append(f"{node(u)} {node(v)} {m} {kind}{k} {fb}")
        lines += [f"bonds {len(bonds)}"] + bonds
        lines += ["starts 1", f"{node(self.start)} {names[self.start[0]]}.{self.start[1]}"]
        with open(path, "w", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")


def _coordinate(node):
    segment, side = node
    return segment.start if side == "h" else segment.end


def _genomic_kind(node):
    """L or R: the genomic side of a node; a segment read right to left (start > end) has its
    start node ('h') on the right."""
    segment, side = node
    return "L" if (side == "h") == (segment.start <= segment.end) else "R"


def _foldback(edge, span):
    """Whether a junction joins two ends of the same genomic kind on one chromosome within span bp."""
    u = next(iter(edge))
    v = _other(edge, u)
    if not (hasattr(u[0], "chr") and hasattr(v[0], "chr")):
        return False
    return is_foldback((u[0].chr, _coordinate(u), _genomic_kind(u)),
                       (v[0].chr, _coordinate(v), _genomic_kind(v)), span)
