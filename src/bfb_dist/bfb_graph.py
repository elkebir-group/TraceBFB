"""BFB graph storage, validation, and output formatting."""

import math

from bfb_dist.history import (
    Operation,
    State,
    end_state,
    start_state,
    validate_sources,
)
from bfb_dist.sequence import Sequence
from bfb_dist.utils import format_operation_line


class BFBGraph:
    """A directed graph containing minimum-length BFB histories.

    Nodes are :class:`~bfb_dist.history.State` objects. ``adjacency`` maps
    each non-end node to a mapping from child states to the
    :class:`~bfb_dist.history.Operation` producing them. Every history
    starts at ``start`` (the empty sequence with no source appended) and
    ends at an *end state*, a state whose sequence is the target. Sources
    are appended at most once, so histories may end with different sets of
    appended sources; ``end`` is the end state that appends every source.
    The graph and its cached results must not be mutated.
    """

    def __init__(self, adjacency, sources, target, bfb_distance=None):
        """Initialize a graph with mandatory adjacency and endpoint data."""
        self.adjacency = adjacency
        self.sources = validate_sources(sources)
        self.target = target
        self.start = start_state()
        self.end = end_state(self.sources, target)
        if bfb_distance is not None:
            self.bfb_distance = bfb_distance

    def is_end(self, state):
        """Return whether ``state`` is an end state, i.e. spells the target."""
        return state.sequence == self.target

    @property
    def ends(self):
        """Return the set of end states of the graph."""
        if len(self.target) == 0:
            return {self.start}
        return {
            child
            for children in self.adjacency.values()
            for child in children
            if self.is_end(child)
        }

    def calculate_bfb_distance(self):
        """Calculate, cache, and return the length of a shortest history.

        Returns ``math.inf`` when no end state is reachable. Valid BFB
        graphs have equal-length start-to-end paths, so following any child
        is sufficient.
        """
        if hasattr(self, "bfb_distance"):
            return self.bfb_distance

        node = self.start
        distance = 0
        while not self.is_end(node):
            children = self.adjacency.get(node)
            if not children:
                distance = math.inf
                break
            node = next(iter(children))
            distance += 1
        self.bfb_distance = distance
        return self.bfb_distance

    @classmethod
    def from_text(cls, input_file):
        """Read and validate a graph from an iterable text stream.

        Raises:
            ValueError: If the graph syntax, operations, or structure are
                invalid.
        """
        lines = iter(input_file)

        sequences = {}

        def parse_sequence(text):
            """Parse and intern one serialized sequence."""
            text = text.strip()
            if text not in sequences:
                sequences[text] = Sequence.from_string(text)
            return sequences[text]

        def parse_state(text):
            """Parse one serialized state."""
            return State.from_string(text, parse_sequence)

        line = next(lines, None)
        if line is None:
            raise ValueError("Invalid graph format: missing target")
        try:
            target = parse_sequence(line.rstrip("\r\n"))
        except ValueError as error:
            raise ValueError(
                "Invalid graph format on line 1: invalid target"
            ) from error

        line = next(lines, None)
        if line is None:
            raise ValueError("Invalid graph format: missing source line")
        line = line.rstrip("\r\n")
        try:
            sources = tuple(
                parse_sequence(field) for field in line.split("\t") if field
            )
        except ValueError as error:
            raise ValueError(
                "Invalid graph format on line 2: invalid source"
            ) from error
        if any(len(source) == 0 for source in sources):
            raise ValueError("Invalid graph format on line 2: empty source")

        adjacency = {}

        for line_number, line in enumerate(lines, start=3):
            line = line.rstrip("\r\n")
            if not line:
                continue

            fields = line.split("\t")
            try:
                node = parse_state(fields[0])
            except ValueError as error:
                raise ValueError(
                    f"Invalid graph format on line {line_number}: "
                    "invalid node"
                ) from error
            if node in adjacency:
                raise ValueError(
                    "Invalid graph format: multiple rows have the same "
                    f"node {node}"
                )

            children = {}
            for serialized_edge in fields[1:]:
                try:
                    child, operation = cls._parse_edge(
                        serialized_edge, parse_state
                    )
                except ValueError as error:
                    raise ValueError(
                        f"Invalid graph format on line {line_number}: "
                        "an edge must be child,kind,x,y or child,t,source"
                    ) from error
                if child in children:
                    raise ValueError(
                        f"Invalid graph format on line {line_number}: "
                        f"multiple edges connect {node} to {child}"
                    )
                children[child] = operation

            adjacency[node] = children

        graph = cls(adjacency, sources, target)
        graph._validate_loaded_graph()
        return graph

    @staticmethod
    def _parse_edge(serialized_edge, parse_state):
        """Parse ``child,kind,x,y`` or ``child,t,source``."""
        fields = serialized_edge.rsplit(",", 3)
        if len(fields) == 4 and fields[1] in ("p", "i", "f"):
            child_text, operation_fields = fields[0], fields[1:]
        else:
            fields = serialized_edge.rsplit(",", 2)
            if len(fields) != 3 or fields[1] != "t":
                raise ValueError(f"invalid edge {serialized_edge!r}")
            child_text, operation_fields = fields[0], fields[1:]
        return parse_state(child_text), Operation.from_fields(operation_fields)

    def write_text(self, output):
        """Write the graph in BFB-dist adjacency-list format."""
        output.write(f"{self.target}\n")
        output.write("\t".join(map(str, self.sources)) + "\n")
        for node in sorted(self.adjacency, key=str):
            fields = [str(node)]
            for child, operation in sorted(
                self.adjacency[node].items(),
                key=lambda item: str(item[0]),
            ):
                fields.append(f"{child},{operation}")
            output.write("\t".join(fields) + "\n")

    def get_paths(self):
        """Enumerate, cache, and return all start-to-end state paths."""
        if hasattr(self, "paths"):
            return self.paths

        if self.is_end(self.start):
            self.paths = [[self.start]]
            return self.paths
        if not self.adjacency.get(self.start):
            self.paths = []
            return self.paths

        self.paths = []
        ordered_children = {
            node: sorted(children, key=str)
            for node, children in self.adjacency.items()
        }
        path = [self.start]

        def visit(node):
            """Append every path reachable from ``node`` to the cache."""
            if self.is_end(node):
                self.paths.append(path.copy())
                return
            for child in ordered_children[node]:
                path.append(child)
                visit(child)
                path.pop()

        visit(self.start)
        return self.paths

    def write_paths(self, colored=True):
        """Return all histories, one sequence per nonempty line.

        The empty start sequence is omitted, so a history of length ``k``
        occupies ``k`` lines: the sequences produced by its first ``k - 1``
        cycles, each colored by the cycle applied to it when ``colored`` is
        true, followed by the target. Alternative histories are separated
        by a blank line.
        """
        formatted_paths = []
        for path in self.get_paths():
            lines = []
            for index in range(1, len(path) - 1):
                node, child = path[index], path[index + 1]
                lines.append(
                    format_operation_line(
                        node.sequence,
                        self.adjacency[node][child],
                        colored=colored,
                    )
                )
            if len(path) > 1:
                lines.append(str(path[-1].sequence))
            if lines:
                formatted_paths.append("\n".join(lines))

        return (
            "\n\n".join(formatted_paths) + "\n"
            if formatted_paths
            else ""
        )

    def write_dot(self, output):
        """Write a Graphviz DOT representation with cycle-kind edge labels."""
        output.write("digraph BFB_dist {\n")
        nodes = sorted(self.adjacency, key=str)
        for node in nodes:
            output.write(f'    "{node}";\n')
        for node in nodes:
            for child, operation in sorted(
                self.adjacency[node].items(),
                key=lambda item: str(item[0]),
            ):
                output.write(
                    f'    "{node}" -> "{child}" [label="{operation.kind}"];\n'
                )
        output.write("}\n")

    def _validate_loaded_graph(self):
        """Validate graph data read from an untrusted text file.

        Validation covers endpoints, declared nodes, BFB operations,
        reachability, terminal behavior, and equal path lengths.
        """
        if self.start not in self.adjacency:
            raise ValueError(
                "Invalid graph format: start state has no adjacency row"
            )
        if self.is_end(self.start):
            if len(self.adjacency) != 1 or self.adjacency[self.start]:
                raise ValueError(
                    "Invalid graph format: a zero-distance graph must "
                    "contain only the start state without edges"
                )
            return

        if not self.adjacency[self.start]:
            if len(self.adjacency) != 1:
                raise ValueError(
                    "Invalid graph format: an unreachable graph must "
                    "contain only the start state without edges"
                )
            return

        for node in self.adjacency:
            if self.is_end(node):
                raise ValueError(
                    "Invalid graph format: end state cannot have an "
                    f"adjacency row: {node}"
                )

        # Valid BFB edges strictly increase sequence length, so the graph
        # cannot contain cycles.
        reachable = {self.start}
        depths = {self.start: 0}
        pending = [self.start]
        while pending:
            node = pending.pop()
            children = self.adjacency[node]
            if not children:
                raise ValueError(
                    "Invalid graph format: non-end node has no "
                    f"outgoing edge: {node}"
                )
            for child, operation in children.items():
                if not self.is_end(child) and child not in self.adjacency:
                    raise ValueError(
                        "Invalid graph format: edge points to an "
                        f"undeclared node: {child}"
                    )
                try:
                    produced = node.apply(operation, self.sources)
                except ValueError as error:
                    raise ValueError(
                        "Invalid graph format: edge operation is not "
                        f"applicable: {node} -> {child}"
                    ) from error
                if produced != child:
                    raise ValueError(
                        "Invalid graph format: edge operation does not "
                        f"generate its child: {node} -> {child}"
                    )

                child_depth = depths[node] + 1
                if child in depths and depths[child] != child_depth:
                    raise ValueError(
                        "Invalid graph format: start-to-end paths "
                        "have different lengths"
                    )
                depths[child] = child_depth

                if child not in reachable:
                    reachable.add(child)
                    if not self.is_end(child):
                        pending.append(child)

        nodes = set(self.adjacency) | self.ends
        if reachable != nodes:
            raise ValueError(
                "Invalid graph format: graph contains nodes that are "
                "not reachable from the start state"
            )
