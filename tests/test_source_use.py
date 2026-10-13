"""Sources may be reused, histories minimize t-BFB cycles first, and t-BFB
cycles may keep the matched prefix."""

import math
import random
import unittest

from bfbdist.distance import BFB_dist, _bounded_search
from bfbdist.nbfb import NBFB_dist
from bfbdist.sequence import LetterAsSegment, Sequence
from tests.helpers import assert_valid_generated_graph, letters

ALGORITHMS = (BFB_dist, NBFB_dist)


class SourceReuseTests(unittest.TestCase):
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
                self.assertEqual(graph.translocation_count(), 1)

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
                self.assertEqual(graph.translocation_count(), 3)

    def test_fewest_translocations_beat_shortest_length(self):
        # The shortest history to A b a B A appends both sources (3 cycles),
        # but a history with a single t-BFB exists (4 cycles) and wins.
        sources = [letters("AB"), letters("BA")]
        target = letters("AbaBA")
        graph = BFB_dist(sources, target)
        self.assert_valid_generated_graph(graph, 4)
        self.assertEqual(graph.translocation_count(), 1)
        _, start_backward, shortest = _bounded_search(sources, target, len(target))
        self.assertEqual((start_backward.translocations, shortest), (2, 3))

    def test_source_may_be_reused(self):
        # A tandem copy costs one t-BFB per copy; the target has no inverted
        # segment, so no fold-back history exists.
        sources, target = [letters("AB")], letters("ABAB")
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                graph = algorithm(sources, target)
                self.assert_valid_generated_graph(graph, 2)
                self.assertEqual(graph.translocation_count(), 2)

    def test_listing_a_source_twice_changes_nothing(self):
        sources = [letters("ABCD"), letters("X")]
        target = letters("ABCDXcbBCDXcb")
        once = BFB_dist(sources, target)
        twice = BFB_dist(sources + sources, target)
        self.assertEqual(once.calculate_bfb_distance(), twice.calculate_bfb_distance())
        self.assertEqual(once.translocation_count(), twice.translocation_count())


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

    def test_general_distance_is_at_most_nested(self):
        rng = random.Random(20260916)
        alphabet = "ABab"
        for _ in range(150):
            target = letters("".join(rng.choice(alphabet) for _ in range(rng.randint(1, 7))))
            sources = [
                letters("".join(rng.choice(alphabet) for _ in range(rng.randint(1, 3))))
                for _ in range(rng.randint(1, 3))
            ]
            with self.subTest(sources=sources, target=target):
                def label(graph):
                    return (graph.translocation_count(), graph.calculate_bfb_distance())
                slow = label(NBFB_dist(sources, target))
                general = label(BFB_dist(sources, target))
                self.assertLessEqual(general, slow)


if __name__ == "__main__":
    unittest.main()
