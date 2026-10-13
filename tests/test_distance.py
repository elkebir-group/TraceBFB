from collections import defaultdict
import math
import unittest

from bfbdist.distance import BFB_dist, _relabel_forward
from bfbdist.history import Operation, State, start_state
from tests.helpers import assert_valid_generated_graph, letters


class DistanceTests(unittest.TestCase):
    def assert_valid_generated_graph(self, graph, expected_distance):
        assert_valid_generated_graph(self, graph, expected_distance)

    def test_single_source_counts_the_opening_translocation(self):
        source = letters("AB")

        graph = BFB_dist([source], source)

        self.assert_valid_generated_graph(graph, 1)
        self.assertEqual(
            graph.adjacency[graph.start][State(graph.target, 1)],
            Operation.translocation(0),
        )

    def test_one_fold_back_after_the_source(self):
        source = letters("AB")
        target = source.bfb_target(0, 2, False)

        graph = BFB_dist([source], target)

        self.assert_valid_generated_graph(graph, 2)

    def test_zero_distance_requires_empty_target_and_no_source(self):
        graph = BFB_dist([], letters(""))

        self.assert_valid_generated_graph(graph, 0)
        self.assertEqual(graph.get_paths(), [[start_state()]])

    def test_unreachable_targets(self):
        cases = [
            ([letters("AB")], letters("C")),
            ([letters("ABC")], letters("AB")),
        ]
        for sources, target in cases:
            with self.subTest(sources=sources, target=target):
                graph = BFB_dist(sources, target)
                self.assert_valid_generated_graph(graph, math.inf)
                self.assertEqual(graph.get_paths(), [])

    def test_empty_source_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "empty"):
            BFB_dist([letters("")], letters("A"))

    def test_empty_r_yields_a_p_bfb_operation(self):
        source = letters("AB")
        target = source.bfb_target(0, len(source), True)

        graph = BFB_dist([source], target)

        self.assert_valid_generated_graph(graph, 2)
        self.assertEqual(
            graph.adjacency[State(source, 1)][State(graph.target, 1)],
            Operation("p", 0, 2),
        )

    def test_nonempty_r_yields_i_and_f_bfb_operations(self):
        source = letters("AB")
        for is_forward_deletion, kind in ((False, "i"), (True, "f")):
            with self.subTest(kind=kind):
                target = source.bfb_target(0, 1, is_forward_deletion)
                graph = BFB_dist([source], target)
                self.assert_valid_generated_graph(graph, 2)
                self.assertEqual(
                    graph.adjacency[State(source, 1)][State(graph.target, 1)],
                    Operation(kind, 0, 1),
                )

    def test_all_minimum_paths_are_retained(self):
        graph = BFB_dist([letters("AB")], letters("AbaAbaA"))

        self.assert_valid_generated_graph(graph, 4)
        self.assertEqual(len(graph.get_paths()), 2)

    def test_longer_paths_are_not_included(self):
        source = letters("AB")
        target = source.bfb_target(0, 2, False)

        graph = BFB_dist([source], target)

        self.assert_valid_generated_graph(graph, 2)
        self.assertEqual(
            graph.get_paths(),
            [[graph.start, State(source, 1), State(graph.target, 1)]],
        )

    def test_paper_examples_are_one_more_than_the_fold_back_count(self):
        # Leftover-separation and prefix-separation examples of the paper,
        # whose single-source fold-back counts are 3 and 4.
        self.assert_valid_generated_graph(
            BFB_dist([letters("AAba")], letters("AAAAAABaa")), 4
        )
        self.assert_valid_generated_graph(
            BFB_dist([letters("AAB")], letters("AABbAbAAB")), 5
        )

    def test_two_sources_in_either_order(self):
        target = letters("ABbC")

        first = BFB_dist([letters("AB"), letters("C")], target)
        second = BFB_dist([letters("C"), letters("AB")], target)

        self.assert_valid_generated_graph(first, 3)
        self.assert_valid_generated_graph(second, 3)
        self.assertEqual(
            [op.kind for op in self._path_operations(first)], ["t", "p", "t"]
        )
        self.assertEqual(
            [op.source for op in self._path_operations(first)], [0, -1, 1]
        )
        self.assertEqual(
            [op.source for op in self._path_operations(second)], [1, -1, 0]
        )

    def test_equal_sources_are_distinguished_by_index(self):
        graph = BFB_dist([letters("A"), letters("A")], letters("AaA"))

        self.assert_valid_generated_graph(graph, 3)
        # Fewest t-BFB cycles first: one copy followed by two p-BFB cycles
        # (A -> A a -> A a A) beats appending both copies. Equal sources
        # are interchangeable, so the parallel t-BFB edges collapse to one
        # path that appends the first copy.
        self.assertEqual(graph.translocation_count(), 1)
        self.assertEqual(len(graph.get_paths()), 1)
        self.assertEqual(self._path_operations(graph)[0], Operation.translocation(0))

    def test_source_may_be_appended_after_fold_backs(self):
        graph = BFB_dist([letters("AB"), letters("C")], letters("ABCcbb"))

        self.assert_valid_generated_graph(graph, 4)

    def test_relabeling_converts_reversed_counts_and_drops_stray_nodes(self):
        # Backward labels count reversed t-BFB cycles: 0 at the target, 1 at
        # the empty sequence. Relabeling turns them into forward counts.
        sources = (letters("AB"),)
        target = State(letters("ABbA"), 0)
        intermediate = State(letters("ABb"), 0)
        source = State(letters("AB"), 0)
        start_backward = State(letters(""), 1)
        stray = State(letters("Z"), 0)
        G = defaultdict(dict)
        G[start_backward][source] = Operation.translocation(0)
        G[source][intermediate] = Operation("i", 0, 1)
        G[intermediate][target] = Operation("f", 1, 2)
        G[stray][State(letters("Zz"), 0)] = Operation("p", 0, 1)

        adjacency = _relabel_forward(G, start_backward, 1)

        forward_source = State(letters("AB"), 1)
        self.assertEqual(
            set(adjacency),
            {start_state(), forward_source, State(letters("ABb"), 1)},
        )
        self.assertEqual(
            adjacency[start_state()],
            {forward_source: Operation.translocation(0)},
        )
        self.assertEqual(
            adjacency[State(letters("ABb"), 1)],
            {State(letters("ABbA"), 1): Operation("f", 1, 2)},
        )
        self.assertNotIn(stray, adjacency)

    def _path_operations(self, graph):
        path = graph.get_paths()[0]
        return [
            graph.adjacency[node][child] for node, child in zip(path, path[1:])
        ]


if __name__ == "__main__":
    unittest.main()
