import math
import random
import unittest

from bfb_dist.distance import BFB_dist
from bfb_dist.history import start_state
from bfb_dist.nbfb import LCP_monotone_dist, NBFB_dist
from bfb_dist.nbfb_fast import NBFB_fast_dist
from bfb_dist.simulation import _set_seed, simulate_history
from tests.helpers import assert_valid_generated_graph, letters


ALGORITHMS = (NBFB_dist, NBFB_fast_dist, LCP_monotone_dist)


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
            ([letters("A")], letters("")),
            ([letters("AB"), letters("C")], letters("ABb")),
            ([letters("C"), letters("AB")], letters("CbAB")),
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

        fast = NBFB_fast_dist([letters("AB")], letters("AbaAbaA"))
        self.assert_valid_generated_graph(fast, 4)
        self.assertEqual(len(fast.get_paths()), 1)

    def test_leftover_separation_example(self):
        # s = a a bar-b bar-a, t = a^6 b bar-a bar-a: LCP-monotone histories
        # exist, nested ones do not.
        sources, target = [letters("AAba")], letters("AAAAAABaa")

        self.assert_valid_generated_graph(LCP_monotone_dist(sources, target), 4)
        self.assert_valid_generated_graph(NBFB_dist(sources, target), math.inf)
        self.assert_valid_generated_graph(
            NBFB_fast_dist(sources, target), math.inf
        )

    def test_prefix_separation_example(self):
        sources, target = [letters("AAB")], letters("AABbAbAAB")
        self.assertTrue(sources[0].is_prefix_of(target))

        self.assert_valid_generated_graph(LCP_monotone_dist(sources, target), 5)
        self.assert_valid_generated_graph(NBFB_dist(sources, target), math.inf)
        self.assert_valid_generated_graph(
            NBFB_fast_dist(sources, target), math.inf
        )

    def test_translocation_after_fold_backs(self):
        sources, target = [letters("AB"), letters("C")], letters("ABbC")

        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                self.assert_valid_generated_graph(algorithm(sources, target), 3)

    def test_source_strings_are_exempt_from_foldability(self):
        # The source A B C has tail B C with respect to the target A c A B C,
        # and c b does not occur in the target. The history
        # t-BFB, f-BFB (p = eps, q = A, r = B C), f-BFB (p = A, q = c, r = b a)
        # is nested only because strings produced by t-BFB cycles are exempt.
        sources, target = [letters("ABC")], letters("AcABC")

        for algorithm in ALGORITHMS:
            with self.subTest(algorithm=algorithm.__name__):
                self.assert_valid_generated_graph(algorithm(sources, target), 3)

    def test_fast_agrees_with_search_on_short_deletion_histories(self):
        # Histories whose f-BFB overhangs are shorter than the preceding
        # cycle's gain are nested (paper, Short deletions), so NBFB is at
        # most the simulated length; NBFB_fast must agree with NBFB and
        # BFB can only be smaller.
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
                nested = distance(NBFB_dist, sources, target)
                fast = distance(NBFB_fast_dist, sources, target)
                general = distance(BFB_dist, sources, target)
                self.assertEqual(nested, fast)
                self.assertLessEqual(nested, cycles)
                self.assertLessEqual(general, nested)
        self.assertGreater(checked, 50)

    def test_fast_agrees_with_search_on_unrestricted_histories(self):
        rng = random.Random(11)
        for trial in range(80):
            n_sources = rng.randint(1, 2)
            sources = [
                letters("".join(rng.choice("AB") for _ in range(rng.randint(1, 3))))
                for _ in range(n_sources)
            ]
            cycles = rng.randint(n_sources, n_sources + 3)
            _set_seed(2000 + trial)
            steps, _ = simulate_history(sources, cycles)
            target = steps[-1][0].sequence
            if len(target) > 14:
                continue
            with self.subTest(sources=sources, target=target):
                self.assertEqual(
                    distance(NBFB_dist, sources, target),
                    distance(NBFB_fast_dist, sources, target),
                )


if __name__ == "__main__":
    unittest.main()
