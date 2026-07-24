"""BFB graph storage, validation, and output formatting."""

import math

from bfb_dist.sequence import Sequence
from bfb_dist.utils import format_sequence


class BFBGraph:
    """A directed graph containing minimum-length BFB paths.

    ``adjacency`` maps each non-target node to a mapping from child sequences
    to ``(x, y, is_forward_deletion)`` operations. The graph and its cached
    results must not be mutated.
    """

    def __init__(self, adjacency, source, target, bfb_distance=None):
        """Initialize a graph with mandatory adjacency and endpoint data."""
        self.adjacency = adjacency
        self.source, self.target = source, target
        if bfb_distance is not None:
            self.bfb_distance = bfb_distance

    def calculate_bfb_distance(self):
        """Calculate, cache, and return the source-to-target distance.

        Returns ``math.inf`` when the target is unreachable. Valid BFB graphs
        have equal-length source-to-target paths, so following any child is
        sufficient.
        """
        if hasattr(self, "bfb_distance"):
            return self.bfb_distance

        node = self.source
        distance = 0
        while node != self.target:
            children = self.adjacency.get(node)
            if not children:
                distance = math.inf
                break
            # Every source-to-target path has the same length, so any child
            # can be retrieved from the adjacency dictionary.
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
        line = next(lines, None)
        if line is None:
            raise ValueError(
                "Invalid graph format: missing source and target"
            )

        line = line.rstrip("\r\n")
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(
                "Invalid graph format on line 1: "
                "expected source and target"
            )

        sequences = {}

        def parse_sequence(text):
            """Parse and intern one serialized sequence."""
            text = text.strip()
            if text not in sequences:
                sequences[text] = Sequence.from_string(text)
            return sequences[text]

        try:
            source, target = map(parse_sequence, fields)
        except ValueError as error:
            raise ValueError(
                "Invalid graph format on line 1: "
                "invalid source or target"
            ) from error

        adjacency = {}

        for line_number, line in enumerate(lines, start=2):
            line = line.rstrip("\r\n")
            if not line:
                continue

            fields = line.split("\t")
            try:
                node = parse_sequence(fields[0])
            except ValueError as error:
                raise ValueError(
                    f"Invalid graph format on line {line_number}: "
                    "invalid source node"
                ) from error
            if node in adjacency:
                raise ValueError(
                    "Invalid graph format: multiple rows have the same "
                    f"source node {node}"
                )

            children = {}
            for serialized_edge in fields[1:]:
                try:
                    child_text, x_text, y_text, forward_text = (
                        serialized_edge.rsplit(",", 3)
                    )
                    child = parse_sequence(child_text)
                    x, y = int(x_text), int(y_text)
                except ValueError as error:
                    raise ValueError(
                        f"Invalid graph format on line {line_number}: "
                        "an edge must be target,x,y,is_forward_deletion"
                    ) from error
                if forward_text not in {"0", "1"}:
                    raise ValueError(
                        f"Invalid graph format on line {line_number}: "
                        "is_forward_deletion must be 0 or 1"
                    )
                if not 0 <= x < y <= len(node):
                    raise ValueError(
                        f"Invalid graph format on line {line_number}: "
                        f"edge coordinates ({x}, {y}) are out of bounds"
                    )
                if child in children:
                    raise ValueError(
                        f"Invalid graph format on line {line_number}: "
                        f"multiple edges connect {node} to {child}"
                    )
                children[child] = (x, y, forward_text == "1")

            adjacency[node] = children

        graph = cls(adjacency, source, target)
        graph._validate_loaded_graph()
        return graph

    def write_text(self, output):
        """Write the graph in BFB-dist adjacency-list format."""
        output.write(f"{self.source}\t{self.target}\n")
        for node in sorted(self.adjacency, key=str):
            fields = [str(node)]
            for child, operation in sorted(
                self.adjacency[node].items(),
                key=lambda item: str(item[0]),
            ):
                x, y, is_forward_deletion = operation
                fields.append(
                    f"{child},{x},{y},{int(is_forward_deletion)}"
                )
            output.write("\t".join(fields) + "\n")

    def get_paths(self):
        """Enumerate, cache, and return all source-to-target node paths."""
        if hasattr(self, "paths"):
            return self.paths

        if self.source == self.target:
            self.paths = [[self.source]]
            return self.paths
        if not self.adjacency.get(self.source):
            self.paths = []
            return self.paths

        self.paths = []
        ordered_children = {
            node: sorted(children, key=str)
            for node, children in self.adjacency.items()
        }
        path = [self.source]

        def visit(node):
            """Append every path reachable from ``node`` to the cache."""
            if node == self.target:
                self.paths.append(path.copy())
                return
            for child in ordered_children[node]:
                path.append(child)
                visit(child)
                path.pop()

        visit(self.source)
        return self.paths

    def write_paths(self, colored=True):
        """Return all paths as complete sequences, one per nonempty line.

        Alternative paths are separated by a blank line. When ``colored`` is
        true, source and intermediate sequences contain ANSI-colored P, Q, and
        R regions.
        """
        formatted_paths = []
        for path in self.get_paths():
            lines = []
            for index in range(len(path) - 1):
                node, child = path[index], path[index + 1]
                if colored:
                    x, y, _ = self.adjacency[node][child]
                    lines.append(
                        format_sequence(
                            node.decompose(x, y),
                            colored=True,
                        )
                    )
                else:
                    lines.append(str(node))
            lines.append(str(path[-1]))
            formatted_paths.append("\n".join(lines))

        return (
            "\n\n".join(formatted_paths) + "\n"
            if formatted_paths
            else ""
        )

    def write_dot(self, output):
        """Write a Graphviz DOT representation without edge labels."""
        output.write("digraph BFB_dist {\n")
        nodes = sorted(self.adjacency, key=str)
        for node in nodes:
            output.write(f'    "{node}";\n')
        for node in nodes:
            for child in sorted(self.adjacency[node], key=str):
                output.write(f'    "{node}" -> "{child}";\n')
        output.write("}\n")

    def _validate_loaded_graph(self):
        """Validate graph data read from an untrusted text file.

        Validation covers endpoints, declared nodes, BFB operations,
        reachability, terminal behavior, and equal path lengths.
        """
        # The source is always the only node in an edge-free graph.
        if self.source not in self.adjacency:
            raise ValueError(
                "Invalid graph format: source has no adjacency row"
            )
        if self.source == self.target:
            if (
                len(self.adjacency) != 1
                or self.adjacency[self.source]
            ):
                raise ValueError(
                    "Invalid graph format: a zero-distance graph must "
                    "contain only the source without edges"
                )
            return

        # An unreachable target is represented by the source alone.
        if not self.adjacency[self.source]:
            if len(self.adjacency) != 1:
                raise ValueError(
                    "Invalid graph format: an unreachable graph must "
                    "contain only the source without edges"
                )
            return

        # A reachable target has no adjacency row because it has no children.
        if self.target in self.adjacency:
            raise ValueError(
                "Invalid graph format: target cannot have an adjacency row"
            )

        # Traverse the reachable graph while validating its nodes and edges.
        # Valid BFB edges strictly increase sequence length, so the resulting
        # graph cannot contain cycles.
        reachable = {self.source}
        depths = {self.source: 0}
        pending = [self.source]
        while pending:
            node = pending.pop()
            children = self.adjacency[node]
            if not children:
                raise ValueError(
                    "Invalid graph format: non-target node has no "
                    f"outgoing edge: {node}"
                )
            for child, operation in children.items():
                if child != self.target and child not in self.adjacency:
                    raise ValueError(
                        "Invalid graph format: edge points to an "
                        f"undeclared node: {child}"
                    )
                if node.bfb_target(*operation) != child:
                    raise ValueError(
                        "Invalid graph format: edge operation does not "
                        f"generate its target: {node} -> {child}"
                    )

                child_depth = depths[node] + 1
                if child in depths and depths[child] != child_depth:
                    raise ValueError(
                        "Invalid graph format: source-to-target paths "
                        "have different lengths"
                    )
                depths[child] = child_depth

                if child not in reachable:
                    reachable.add(child)
                    if child != self.target:
                        pending.append(child)

        # Every declared node, including the target, must have been reached.
        # Together with acyclicity and outgoing edges above, this guarantees
        # that every directed path ends at the target.
        nodes = set(self.adjacency) | {self.target}
        if reachable != nodes:
            raise ValueError(
                "Invalid graph format: graph contains nodes that are "
                "not reachable from the source"
            )
