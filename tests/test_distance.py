from collections import defaultdict
import math
import unittest

from bfb_dist.bfb_graph import BFBGraph
from bfb_dist.distance import BFB_dist, _truncate_G
from tests.helpers import letters


class DistanceTests(unittest.TestCase):
    def assert_valid_generated_graph(self, graph, expected_distance):
        self.assertIsInstance(graph, BFBGraph)
        self.assertEqual(graph.calculate_bfb_distance(), expected_distance)

        if expected_distance == 0 or math.isinf(expected_distance):
            self.assertEqual(graph.adjacency, {graph.source: {}})
            return

        self.assertNotIn(graph.target, graph.adjacency)
        for source, children in graph.adjacency.items():
            self.assertTrue(children)
            for target, operation in children.items():
                self.assertEqual(
                    source.bfb_target(*operation),
                    target,
                )

        paths = graph.get_paths()
        self.assertTrue(paths)
        for path in paths:
            self.assertEqual(path[0], graph.source)
            self.assertEqual(path[-1], graph.target)
            self.assertEqual(len(path) - 1, expected_distance)

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
