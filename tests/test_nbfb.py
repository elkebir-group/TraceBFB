import math
import unittest

from bfb_dist.nbfb import NBFB_dist, NBFB_star_dist
from tests.helpers import assert_valid_generated_graph, letters


class NbfbTests(unittest.TestCase):
    def assert_valid_generated_graph(self, graph, expected_distance):
        assert_valid_generated_graph(self, graph, expected_distance)

    def test_zero_distance(self):
        source = letters("AB")

        graph = NBFB_dist(source, source)

        self.assert_valid_generated_graph(graph, 0)
        self.assertEqual(graph.get_paths(), [[source]])

    def test_zero_distance_star(self):
        source = letters("AB")

        graph = NBFB_star_dist(source, source)

        self.assert_valid_generated_graph(graph, 0)
        self.assertEqual(graph.get_paths(), [[source]])

    def test_unreachable_target(self):
        graph = NBFB_dist(letters("AB"), letters("C"))

        self.assert_valid_generated_graph(graph, math.inf)
        self.assertEqual(graph.get_paths(), [])

    def test_unreachable_target_star(self):
        graph = NBFB_star_dist(letters("AB"), letters("C"))

        self.assert_valid_generated_graph(graph, math.inf)
        self.assertEqual(graph.get_paths(), [])

    def test_all_minimum_paths_are_retained(self):
        source = letters("AB")
        target = letters("AbaAbaA")

        graph = NBFB_dist(source, target)

        self.assert_valid_generated_graph(graph, 3)

    def test_leftover_separation_example(self):
        # paper.tex:761-778 (Example leftover-separation).
        # s = a a bar-b bar-a, t = a^6 b bar-a bar-a.
        source = letters("AAba")
        target = letters("AAAAAABaa")

        graph = NBFB_dist(source, target)
        self.assert_valid_generated_graph(graph, 3)

        star_graph = NBFB_star_dist(source, target)
        self.assert_valid_generated_graph(star_graph, math.inf)

    def test_prefix_separation_example(self):
        # paper.tex:790-799, exact value confirmed at the meknote on
        # paper.tex:851 ("NBFB(s,t) = 4 exactly").
        # s = a a b, t = a a b bar-b a bar-b a a b.
        source = letters("AAB")
        target = letters("AABbAbAAB")
        self.assertTrue(source.is_prefix_of(target))

        graph = NBFB_dist(source, target)
        self.assert_valid_generated_graph(graph, 4)

        star_graph = NBFB_star_dist(source, target)
        self.assert_valid_generated_graph(star_graph, math.inf)

    def test_nbfb_star_upper_bounds_nbfb(self):
        # Proposition star-vs-nbfb (paper.tex:750-756): NBFB <= NBFB*
        # always. Small, hand-picked instances -- NBFB's search space is
        # worst-case exponential (paper.tex:481-488), so no randomization.
        cases = [
            (letters("AB"), letters("AB").bfb_target(0, 1, False)),
            (letters("AB"), letters("AbaAbaA")),
            (letters("AAB"), letters("AABbAbAAB")),
            (letters("AAba"), letters("AAAAAABaa")),
        ]
        for source, target in cases:
            with self.subTest(source=source, target=target):
                nbfb = NBFB_dist(source, target).calculate_bfb_distance()
                nbfb_star = NBFB_star_dist(
                    source, target
                ).calculate_bfb_distance()
                self.assertLessEqual(nbfb, nbfb_star)


if __name__ == "__main__":
    unittest.main()
