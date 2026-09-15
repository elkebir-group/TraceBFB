from io import StringIO
import math
import re
import unittest

from bfb_dist.bfb_graph import BFBGraph
from bfb_dist.history import Operation, State, start_state
from tests.helpers import genomic_sequence


ANSI = re.compile(r"\033\[[0-9;]*m")


def linear_graph():
    """Return a three-cycle history t-BFB, i-BFB, f-BFB as a graph."""
    source = genomic_sequence(2)
    sources = (source,)
    start = start_state()
    appended = State(source, 1)
    first_operation = Operation("i", 0, 1)
    intermediate = appended.apply(first_operation, sources)
    second_operation = Operation("f", 1, 2)
    end = intermediate.apply(second_operation, sources)
    adjacency = {
        start: {appended: Operation.translocation(0)},
        appended: {intermediate: first_operation},
        intermediate: {end: second_operation},
    }
    graph = BFBGraph(adjacency, sources, end.sequence)
    return graph, appended, intermediate


class BFBGraphTests(unittest.TestCase):
    def test_constructor_retains_adjacency_and_optional_distance(self):
        graph, _, _ = linear_graph()
        adjacency = graph.adjacency

        self.assertIs(graph.adjacency, adjacency)
        self.assertFalse(hasattr(graph, "bfb_distance"))
        self.assertEqual(graph.start, start_state())
        self.assertEqual(graph.end, State(graph.target, 1))

        cached = BFBGraph(adjacency, graph.sources, graph.target, 9)
        self.assertEqual(cached.calculate_bfb_distance(), 9)

    def test_calculates_and_caches_distance(self):
        graph, _, _ = linear_graph()

        self.assertEqual(graph.calculate_bfb_distance(), 3)
        self.assertEqual(graph.bfb_distance, 3)
        self.assertEqual(graph.calculate_bfb_distance(), 3)

    def test_get_paths_and_formatting(self):
        graph, appended, intermediate = linear_graph()

        self.assertEqual(
            graph.get_paths(),
            [[graph.start, appended, intermediate, graph.end]],
        )
        self.assertIs(graph.get_paths(), graph.paths)

        plain = graph.write_paths(colored=False)
        colored = graph.write_paths(colored=True)
        # The empty start is omitted: one line per cycle.
        self.assertEqual(
            plain.splitlines(),
            [str(state.sequence) for state in graph.get_paths()[0][1:]],
        )
        self.assertEqual(
            ANSI.sub("", colored).splitlines(),
            plain.splitlines(),
        )
        self.assertIn("\033[91m", colored)
        self.assertIn("\033[92m", colored)
        self.assertIn("\033[95m", colored)

    def test_translocation_lines_are_colored_cyan(self):
        graph, _, _ = linear_graph()
        adjacency = {
            graph.start: {State(graph.sources[0], 1): Operation.translocation(0)}
        }
        chain = BFBGraph(adjacency, graph.sources, graph.sources[0])

        colored = chain.write_paths(colored=True)

        self.assertEqual(colored, f"{graph.sources[0]}\n")
        self.assertNotIn("\033[96m", colored)

        two_step = BFBGraph(
            {
                graph.start: {
                    State(graph.sources[0], 1): Operation.translocation(0)
                },
                State(graph.sources[0], 1): {
                    State(graph.sources[0] + graph.sources[0], 0b11):
                        Operation.translocation(1)
                },
            },
            (graph.sources[0], graph.sources[0]),
            graph.sources[0] + graph.sources[0],
        )
        colored = two_step.write_paths(colored=True)
        self.assertIn("\033[96m", colored.splitlines()[0])

    def test_text_round_trip(self):
        graph, _, _ = linear_graph()
        output = StringIO()
        graph.write_text(output)

        loaded = BFBGraph.from_text(StringIO(output.getvalue()))

        self.assertEqual(loaded.sources, graph.sources)
        self.assertEqual(loaded.target, graph.target)
        self.assertEqual(loaded.adjacency, graph.adjacency)
        self.assertNotIn(loaded.end, loaded.adjacency)
        self.assertEqual(loaded.calculate_bfb_distance(), 3)

    def test_text_round_trip_with_two_sources(self):
        first, second = genomic_sequence(1, "chr1"), genomic_sequence(1, "chr2")
        sources = (first, second)
        start = start_state()
        one = State(first, 0b01)
        both = State(first + second, 0b11)
        graph = BFBGraph(
            {
                start: {one: Operation.translocation(0)},
                one: {both: Operation.translocation(1)},
            },
            sources,
            first + second,
        )
        output = StringIO()
        graph.write_text(output)

        loaded = BFBGraph.from_text(StringIO(output.getvalue()))

        self.assertEqual(loaded.sources, sources)
        self.assertEqual(loaded.adjacency, graph.adjacency)
        self.assertEqual(loaded.calculate_bfb_distance(), 2)

    def test_dot_labels_edges_with_cycle_kind(self):
        graph, appended, intermediate = linear_graph()
        output = StringIO()

        graph.write_dot(output)
        dot = output.getvalue()

        self.assertIn(f'"{graph.start}" -> "{appended}" [label="t"];', dot)
        self.assertIn(f'"{appended}" -> "{intermediate}" [label="i"];', dot)
        self.assertIn(f'"{intermediate}" -> "{graph.end}" [label="f"];', dot)

    def test_zero_distance_graph(self):
        start = start_state()
        graph = BFBGraph({start: {}}, (), start.sequence)

        self.assertEqual(graph.calculate_bfb_distance(), 0)
        self.assertEqual(graph.get_paths(), [[start]])
        self.assertEqual(graph.write_paths(False), "")

    def test_unreachable_graph(self):
        source = genomic_sequence(2)
        start = start_state()
        graph = BFBGraph({start: {}}, (source,), source.reverse_complement())

        self.assertEqual(graph.calculate_bfb_distance(), math.inf)
        self.assertEqual(graph.get_paths(), [])
        self.assertEqual(graph.write_paths(), "")


class LoadedGraphValidationTests(unittest.TestCase):
    def setUp(self):
        self.graph, self.appended, self.intermediate = linear_graph()
        self.source = self.graph.sources[0]
        self.start = self.graph.start
        self.end = self.graph.end
        self.header = f"{self.graph.target}\n{self.source}\n"

    def assert_invalid(self, text, message):
        with self.assertRaisesRegex(ValueError, message):
            BFBGraph.from_text(StringIO(text))

    def test_missing_header(self):
        self.assert_invalid("", "missing target")
        self.assert_invalid(f"{self.graph.target}\n", "missing source line")

    def test_missing_start_row(self):
        self.assert_invalid(self.header, "start state has no adjacency row")

    def test_duplicate_row(self):
        self.assert_invalid(
            self.header + f"{self.start}\n{self.start}\n",
            "multiple rows have the same node",
        )

    def test_invalid_operation(self):
        self.assert_invalid(
            self.header + f"{self.start}\t{self.appended},q,0,1\n",
            "an edge must be",
        )

    def test_inapplicable_operation(self):
        self.assert_invalid(
            f"{self.source}\n{self.source}\n"
            f"{self.start}\t{self.appended},p,0,1\n",
            "not applicable",
        )

    def test_end_cannot_have_adjacency_row(self):
        self.assert_invalid(
            f"{self.source}\n{self.source}\n"
            f"{self.start}\t{self.appended},t,0\n"
            f"{self.appended}\n",
            "end state cannot have an adjacency row",
        )

    def test_edge_operation_must_generate_child(self):
        self.assert_invalid(
            f"{self.source}\n{self.source}\n"
            f"{self.start}\t{self.appended},t,1\n",
            "not applicable",
        )
        wrong = State(self.source.reverse_complement(), 1)
        self.assert_invalid(
            f"{wrong.sequence}\n{self.source}\n"
            f"{self.start}\t{wrong},t,0\n",
            "does not generate its child",
        )

    def test_undeclared_intermediate_is_rejected(self):
        self.assert_invalid(
            self.header + f"{self.start}\t{self.appended},t,0\n",
            "undeclared node",
        )

    def test_paths_must_have_equal_lengths(self):
        # A -> A a -> A a A a takes two p-BFB cycles, while
        # A -> A a -> A a A -> A a A a takes three; every edge is valid.
        sources = (genomic_sequence(1),)
        appended = State(sources[0], 1)
        aa = appended.apply(Operation("p", 0, 1), sources)
        aaaa = aa.apply(Operation("p", 0, 2), sources)
        aaa = aa.apply(Operation("p", 1, 2), sources)
        self.assertEqual(aaa.apply(Operation("p", 2, 3), sources), aaaa)
        graph = BFBGraph(
            {
                self.start: {appended: Operation.translocation(0)},
                appended: {aa: Operation("p", 0, 1)},
                aa: {aaaa: Operation("p", 0, 2), aaa: Operation("p", 1, 2)},
                aaa: {aaaa: Operation("p", 2, 3)},
            },
            sources,
            aaaa.sequence,
        )
        serialized = StringIO()
        graph.write_text(serialized)

        self.assert_invalid(serialized.getvalue(), "different lengths")


if __name__ == "__main__":
    unittest.main()
