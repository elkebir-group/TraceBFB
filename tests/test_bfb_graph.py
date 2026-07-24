from io import StringIO
import math
import re
import unittest

from bfb_dist.bfb_graph import BFBGraph
from tests.helpers import genomic_sequence


ANSI = re.compile(r"\033\[[0-9;]*m")


def linear_graph():
    source = genomic_sequence(2)
    first_operation = (0, 1, False)
    intermediate = source.bfb_target(*first_operation)
    second_operation = (1, 2, True)
    target = intermediate.bfb_target(*second_operation)
    adjacency = {
        source: {intermediate: first_operation},
        intermediate: {target: second_operation},
    }
    return BFBGraph(adjacency, source, target), intermediate


class BFBGraphTests(unittest.TestCase):
    def test_constructor_retains_adjacency_and_optional_distance(self):
        graph, _ = linear_graph()
        adjacency = graph.adjacency

        self.assertIs(graph.adjacency, adjacency)
        self.assertFalse(hasattr(graph, "bfb_distance"))

        cached = BFBGraph(adjacency, graph.source, graph.target, 9)
        self.assertEqual(cached.calculate_bfb_distance(), 9)

    def test_calculates_and_caches_distance(self):
        graph, _ = linear_graph()

        self.assertEqual(graph.calculate_bfb_distance(), 2)
        self.assertEqual(graph.bfb_distance, 2)
        self.assertEqual(graph.calculate_bfb_distance(), 2)

    def test_get_paths_and_formatting(self):
        graph, intermediate = linear_graph()

        self.assertEqual(
            graph.get_paths(),
            [[graph.source, intermediate, graph.target]],
        )
        self.assertIs(graph.get_paths(), graph.paths)

        plain = graph.write_paths(colored=False)
        colored = graph.write_paths(colored=True)
        self.assertEqual(
            plain.splitlines(),
            list(map(str, graph.get_paths()[0])),
        )
        self.assertEqual(
            ANSI.sub("", colored).splitlines(),
            plain.splitlines(),
        )
        self.assertIn("\033[91m", colored)
        self.assertIn("\033[92m", colored)
        self.assertIn("\033[95m", colored)
        self.assertNotIn("R: ", plain)
        self.assertNotIn("-> ", plain)

    def test_text_round_trip(self):
        graph, _ = linear_graph()
        output = StringIO()
        graph.write_text(output)

        loaded = BFBGraph.from_text(StringIO(output.getvalue()))

        self.assertEqual(loaded.source, graph.source)
        self.assertEqual(loaded.target, graph.target)
        self.assertEqual(loaded.adjacency, graph.adjacency)
        self.assertNotIn(loaded.target, loaded.adjacency)
        self.assertEqual(loaded.calculate_bfb_distance(), 2)

    def test_dot_uses_sequence_nodes_without_edge_labels(self):
        graph, intermediate = linear_graph()
        output = StringIO()

        graph.write_dot(output)
        dot = output.getvalue()

        self.assertIn(f'"{graph.source}" -> "{intermediate}";', dot)
        self.assertIn(f'"{intermediate}" -> "{graph.target}";', dot)
        self.assertNotIn("label=", dot)

    def test_zero_distance_graph(self):
        source = genomic_sequence(2)
        graph = BFBGraph({source: {}}, source, source)

        self.assertEqual(graph.calculate_bfb_distance(), 0)
        self.assertEqual(graph.get_paths(), [[source]])
        self.assertEqual(graph.write_paths(False), f"{source}\n")

    def test_unreachable_graph(self):
        source = genomic_sequence(2)
        target = source.reverse_complement()
        graph = BFBGraph({source: {}}, source, target)

        self.assertEqual(graph.calculate_bfb_distance(), math.inf)
        self.assertEqual(graph.get_paths(), [])
        self.assertEqual(graph.write_paths(), "")


class LoadedGraphValidationTests(unittest.TestCase):
    def setUp(self):
        self.graph, self.intermediate = linear_graph()
        self.source = self.graph.source
        self.target = self.graph.target

    def assert_invalid(self, text, message):
        with self.assertRaisesRegex(ValueError, message):
            BFBGraph.from_text(StringIO(text))

    def test_missing_header(self):
        self.assert_invalid("", "missing source and target")

    def test_missing_source_row(self):
        self.assert_invalid(
            f"{self.source}\t{self.target}\n",
            "source has no adjacency row",
        )

    def test_duplicate_source_row(self):
        self.assert_invalid(
            f"{self.source}\t{self.target}\n"
            f"{self.source}\n"
            f"{self.source}\n",
            "multiple rows have the same source node",
        )

    def test_invalid_operation_flag(self):
        self.assert_invalid(
            f"{self.source}\t{self.intermediate}\n"
            f"{self.source}\t{self.intermediate},0,1,2\n",
            "is_forward_deletion must be 0 or 1",
        )

    def test_out_of_bounds_operation(self):
        self.assert_invalid(
            f"{self.source}\t{self.intermediate}\n"
            f"{self.source}\t{self.intermediate},0,0,0\n",
            "out of bounds",
        )

    def test_target_cannot_have_adjacency_row(self):
        direct_target = self.source.bfb_target(0, 1, False)
        self.assert_invalid(
            f"{self.source}\t{direct_target}\n"
            f"{self.source}\t{direct_target},0,1,0\n"
            f"{direct_target}\n",
            "target cannot have an adjacency row",
        )

    def test_edge_operation_must_generate_child(self):
        wrong_target = self.source.bfb_target(1, 2, False)
        self.assert_invalid(
            f"{self.source}\t{wrong_target}\n"
            f"{self.source}\t{wrong_target},0,1,0\n",
            "does not generate its target",
        )

    def test_undeclared_intermediate_is_rejected(self):
        self.assert_invalid(
            f"{self.source}\t{self.target}\n"
            f"{self.source}\t{self.intermediate},0,1,0\n",
            "undeclared node",
        )

    def test_source_to_target_paths_must_have_equal_lengths(self):
        source = genomic_sequence(2)
        intermediate = source.bfb_target(1, 2, False)
        target = source.bfb_target(0, 2, False)
        graph = BFBGraph(
            {
                source: {
                    intermediate: (1, 2, False),
                    target: (0, 2, False),
                },
                intermediate: {target: (0, 1, False)},
            },
            source,
            target,
        )
        serialized = StringIO()
        graph.write_text(serialized)

        self.assert_invalid(
            serialized.getvalue(),
            "paths have different lengths",
        )


if __name__ == "__main__":
    unittest.main()
