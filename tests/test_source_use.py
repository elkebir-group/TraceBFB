"""Sources are appended at most once, and t-BFB cycles may keep the prefix."""

import math
import random
import unittest

from bfb_dist.distance import BFB_dist
from bfb_dist.nbfb import LCP_monotone_dist, NBFB_dist
from bfb_dist.nbfb_fast import NBFB_fast_dist
from bfb_dist.sequence import LetterAsSegment, Sequence
from tests.helpers import assert_valid_generated_graph, letters

ALGORITHMS = (BFB_dist, NBFB_dist, NBFB_fast_dist, LCP_monotone_dist)


class AtMostOnceTests(unittest.TestCase):
    def assert_valid_generated_graph(self, graph, expected_distance):
        assert_valid_generated_graph(self, graph, expected_distance)

    def test_empty_target_needs_no_source(self):
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                graph = algorithm([letters("A")], letters(""))
                self.assert_valid_generated_graph(graph, 0)

    def test_unused_source_is_allowed(self):
        sources, target = [letters("AB"), letters("C")], letters("ABb")
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                graph = algorithm(sources, target)
                self.assert_valid_generated_graph(graph, 2)
                for path in graph.get_paths():
                    self.assertEqual(path[-1].used, 0b01)

    def test_unused_sources_do_not_change_distance(self):
        target = letters("ABbaAB")
        base = BFB_dist([letters("AB")], target).calculate_bfb_distance()
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                graph = algorithm(
                    [letters("AB"), letters("ABC"), letters("D")], target
                )
                self.assertEqual(graph.calculate_bfb_distance(), base)

    def test_hg008_blocks_recover_founder(self):
        # Germline Z A B C D E F G P; the target lacks D and has an extra F.
        sources = [
            letters("ZABCDEFGP"),
            letters("ZABC"),
            letters("F"),
            letters("EFGP"),
        ]
        target = letters("ZABCFEFGPgfefcbaB")
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                graph = algorithm(sources, target)
                self.assert_valid_generated_graph(graph, 5)
                for path in graph.get_paths():
                    self.assertEqual(path[-1].used, 0b1110)

    def test_backward_search_relabels_each_unused_set(self):
        # Two different unused sets reach the same distance.
        sources = [letters("AB"), letters("A"), letters("B")]
        target = letters("ABb")
        graph = BFB_dist(sources, target)
        self.assert_valid_generated_graph(graph, 2)
        self.assertEqual({path[-1].used for path in graph.get_paths()}, {0b001})


class RelaxedNestedTests(unittest.TestCase):
    """A t-BFB may keep the matched prefix if the tail stays foldable."""

    def assert_valid_generated_graph(self, graph, expected_distance):
        assert_valid_generated_graph(self, graph, expected_distance)

    def test_apex_insert_at_i_type_fold(self):
        sources, target = [letters("ABCD"), letters("X")], letters("ABCDXcb")
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                self.assert_valid_generated_graph(algorithm(sources, target), 3)

    def test_apex_insert_at_f_type_fold(self):
        # A B C D -> append x -> f-BFB with overhang D x: A B C X d c b.
        sources, target = [letters("ABCD"), letters("x")], letters("ABCXdcb")
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                self.assert_valid_generated_graph(algorithm(sources, target), 3)

    def test_prefix_kept_from_empty_tail(self):
        # A B -> append x -> f-BFB with overhang x: A B X b.
        sources, target = [letters("AB"), letters("x")], letters("ABXb")
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                self.assert_valid_generated_graph(algorithm(sources, target), 3)

    def test_prefix_kept_from_window_tail(self):
        # A B C D -> i-BFB (q = B C, r = D) -> A B C D c b, whose tail b is a
        # window; append y without advancing, then f-BFB with overhang b y:
        # A B C D c Y B C.
        sources = [letters("ABCD"), letters("y")]
        target = letters("ABCDcYBC")
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                graph = algorithm(sources, target)
                self.assert_valid_generated_graph(graph, 4)
                path = graph.get_paths()[0]
                kinds = [
                    graph.adjacency[a][b].kind for a, b in zip(path, path[1:])
                ]
                self.assertEqual(kinds, ["t", "i", "t", "f"])

    def test_fast_and_slow_nested_searches_agree(self):
        rng = random.Random(20260916)
        alphabet = "ABab"
        for _ in range(150):
            target = letters("".join(rng.choice(alphabet) for _ in range(rng.randint(1, 7))))
            sources = [
                letters("".join(rng.choice(alphabet) for _ in range(rng.randint(1, 3))))
                for _ in range(rng.randint(1, 3))
            ]
            with self.subTest(sources=sources, target=target):
                slow = NBFB_dist(sources, target).calculate_bfb_distance()
                fast = NBFB_fast_dist(sources, target).calculate_bfb_distance()
                self.assertEqual(slow, fast)
                general = BFB_dist(sources, target).calculate_bfb_distance()
                self.assertLessEqual(general, slow)


if __name__ == "__main__":
    unittest.main()
