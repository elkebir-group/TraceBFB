import math
import random
import unittest

from bfb_dist.distance import BFB_dist
from bfb_dist.history import start_state
from bfb_dist.nbfb import NBFB_dist, _bounded_search
from bfb_dist.simulation import _set_seed, simulate_history
from tests.helpers import assert_valid_generated_graph, letters


ALGORITHMS = (NBFB_dist,)


def distance(algorithm, sources, target):
    return algorithm(sources, target).calculate_bfb_distance()


class NbfbTests(unittest.TestCase):
    def assert_valid_generated_graph(self, graph, expected_distance):
        assert_valid_generated_graph(self, graph, expected_distance)

    def test_zero_distance(self):
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                graph = algorithm([], letters(""))
                self.assert_valid_generated_graph(graph, 0)
                self.assertEqual(graph.get_paths(), [[start_state()]])

    def test_single_source_equal_to_target(self):
        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                graph = algorithm([letters("AB")], letters("AB"))
                self.assert_valid_generated_graph(graph, 1)

    def test_unreachable_target(self):
        cases = [
            ([letters("AB")], letters("C")),
            ([letters("C"), letters("AB")], letters("CbAB")),
            ([letters("ABC")], letters("AB")),
        ]
        for algorithm in ALGORITHMS:
            for sources, target in cases:
                with self.subTest(
                    algorithm=algorithm.__name__, sources=sources, target=target
                ):
                    graph = algorithm(sources, target)
                    self.assert_valid_generated_graph(graph, math.inf)
                    self.assertEqual(graph.get_paths(), [])

    def test_all_minimum_paths_are_retained(self):
        graph = NBFB_dist([letters("AB")], letters("AbaAbaA"))

        self.assert_valid_generated_graph(graph, 4)
        self.assertEqual(len(graph.get_paths()), 2)

    def test_leftover_example_is_nested(self):
        # s = a a bar-b bar-a, t = a^6 b bar-a bar-a: excluded by the
        # earlier LCP-based definition, nested under the prefix definition
        # (t, i, i, f).
        sources, target = [letters("AAba")], letters("AAAAAABaa")
        graph = NBFB_dist(sources, target)
        self.assert_valid_generated_graph(graph, 4)
        self.assertEqual(graph.translocation_count(), 1)

    def test_prefix_example_needs_one_translocation(self):
        # The source is appended once (t, p, f, f, f); no second t-BFB is
        # needed under the prefix definition.
        sources, target = [letters("AAB")], letters("AABbAbAAB")
        graph = NBFB_dist(sources, target)
        self.assert_valid_generated_graph(graph, 5)
        self.assertEqual(graph.translocation_count(), 1)
        _, _, ends, _ = _bounded_search(sources, target, 1)
        self.assertNotEqual(ends, [])

    def test_translocation_after_fold_backs(self):
        sources, target = [letters("AB"), letters("C")], letters("ABbC")

        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                self.assert_valid_generated_graph(algorithm(sources, target), 3)

    def test_source_suffix_makes_a_string_nested(self):
        # The source A B C has tail B C with respect to the target A c A B C,
        # and c b does not occur in the target. The history
        # t-BFB, f-BFB (p = eps, q = A, r = B C), f-BFB (p = A, q = c, r = b a)
        # is nested because the first string ends with the source.
        sources, target = [letters("ABC")], letters("AcABC")

        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                self.assert_valid_generated_graph(algorithm(sources, target), 3)

    def test_short_deletion_histories_are_nested(self):
        # Histories whose f-BFB overhangs are shorter than the preceding
        # cycle's gain are nested (paper, Short deletions). The simulated
        # history appends every source once, so the optimum, which minimizes
        # t-BFB cycles first and total cycles second, is at most
        # (n_sources, cycles) in that order, and BFB can only be smaller.
        rng = random.Random(3)
        checked = 0
        for trial in range(120):
            n_sources = rng.randint(1, 3)
            sources = [
                letters("".join(rng.choice("ABC") for _ in range(rng.randint(1, 3))))
                for _ in range(n_sources)
            ]
            cycles = rng.randint(n_sources, n_sources + 3)
            _set_seed(1000 + trial)
            steps, _ = simulate_history(sources, cycles, short_deletions=True)
            target = steps[-1][0].sequence
            if len(target) > 14:
                continue
            checked += 1
            with self.subTest(sources=sources, target=target, cycles=cycles):
                def label(algorithm):
                    graph = algorithm(sources, target)
                    return (
                        graph.translocation_count(),
                        graph.calculate_bfb_distance(),
                    )
                nested = label(NBFB_dist)
                general = label(BFB_dist)
                self.assertLessEqual(nested, (n_sources, cycles))
                self.assertLessEqual(general, nested)
        self.assertGreater(checked, 50)


if __name__ == "__main__":
    unittest.main()
