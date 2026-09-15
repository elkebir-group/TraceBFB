from collections import defaultdict
import math
import unittest

from bfb_dist.distance import BFB_dist, _truncate_G
from tests.helpers import assert_valid_generated_graph, letters


class DistanceTests(unittest.TestCase):
    def assert_valid_generated_graph(self, graph, expected_distance):
        assert_valid_generated_graph(self, graph, expected_distance)

    def test_direct_python_api_returns_bfb_graph(self):
        source = letters("AB")
        target = source.bfb_target(0, 1, False)

        graph = BFB_dist(source, target)

        self.assert_valid_generated_graph(graph, 1)

    def test_zero_distance(self):
        source = letters("AB")

        graph = BFB_dist(source, source)

        self.assert_valid_generated_graph(graph, 0)
        self.assertEqual(graph.get_paths(), [[source]])

    def test_unreachable_target(self):
        graph = BFB_dist(letters("AB"), letters("C"))

        self.assert_valid_generated_graph(graph, math.inf)
        self.assertEqual(graph.get_paths(), [])

    def test_empty_r_is_canonicalized_to_not_forward_deletion(self):
        source = letters("AB")
        target = source.bfb_target(0, len(source), True)

        graph = BFB_dist(source, target)

        self.assert_valid_generated_graph(graph, 1)
        self.assertEqual(
            graph.adjacency[source][target],
            (0, len(source), False),
        )

    def test_nonempty_r_retains_forward_deletion(self):
        source = letters("AB")
        target = source.bfb_target(0, 1, True)

        graph = BFB_dist(source, target)

        self.assert_valid_generated_graph(graph, 1)
        self.assertEqual(graph.adjacency[source][target], (0, 1, True))

    def test_all_minimum_paths_are_retained(self):
        source = letters("AB")
        target = letters("AbaAbaA")

        graph = BFB_dist(source, target)

        self.assert_valid_generated_graph(graph, 3)
        self.assertEqual(len(graph.get_paths()), 2)

    def test_longer_paths_are_not_included(self):
        source = letters("AB")
        target = source.bfb_target(0, 2, False)

        graph = BFB_dist(source, target)

        self.assert_valid_generated_graph(graph, 1)
        self.assertEqual(graph.get_paths(), [[source, target]])

    def test_truncation_removes_nodes_not_reachable_from_source(self):
        source = letters("AB")
        intermediate = source.bfb_target(0, 1, False)
        target = intermediate.bfb_target(1, 2, True)
        stray = letters("Z")
        G = defaultdict(dict)
        G[source][intermediate] = (0, 1, False)
        G[intermediate][target] = (1, 2, True)
        G[stray][letters("Zz")] = (0, 1, False)

        adjacency = _truncate_G(G, source)

        self.assertEqual(set(adjacency), {source, intermediate})
        self.assertNotIn(target, adjacency)
        self.assertNotIn(stray, adjacency)


if __name__ == "__main__":
    unittest.main()
