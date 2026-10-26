"""Graph construction (`bfbdist graph`) without BAMs: junction filters, segments,
the fit, background removal, and the graph file."""

import tempfile
import unittest
from collections import Counter
from pathlib import Path

from bfbdist.graph_construction import (
    DEFAULTS, Locus, cut_segments, fit, segment_ids, remove_background, somatic_junctions,
)
from bfbdist.graph_file import GraphFile, GraphSegment

WINDOW = [("chr4", 1000, 9000)]


def call(a, b, n):
    return ((("chr4",) + a), (("chr4",) + b), n)


class JunctionTests(unittest.TestCase):
    def test_filters(self):
        tumor = [
            call((3000, "L"), (3500, "L"), 16),   # fold-back: kept
            call((5000, "R"), (5800, "R"), 17),   # fold-back: kept
            call((4000, "R"), (4300, "L"), 30),   # small deletion: dropped
            call((6000, "R"), (7000, "L"), 3),    # low support: dropped
            call((2000, "R"), (2600, "R"), 9),    # germline: dropped
            call((8000, "R"), (20000, "L"), 9),   # one end outside: dropped
        ]
        normal = [call((2050, "R"), (2620, "R"), 4)]
        js = somatic_junctions(tumor, normal, WINDOW, DEFAULTS)
        self.assertEqual([(j["a"][1], j["b"][1]) for j in js], [(3000, 3500), (5000, 5800)])
        self.assertEqual([j["name"] for j in js], ["j0", "j1"])

    def test_support_is_best_nearby_call(self):
        tumor = [call((3000, "L"), (3500, "L"), 16), call((3510, "L"), (3020, "L"), 40)]
        js = somatic_junctions(tumor, [], WINDOW, DEFAULTS)
        self.assertEqual([j["support"] for j in js], [40, 40])


class SegmentTests(unittest.TestCase):
    def test_cut_and_name(self):
        js = somatic_junctions([call((3000, "L"), (3500, "L"), 16)], [], WINDOW, DEFAULTS)
        segs = cut_segments(WINDOW, js)
        self.assertEqual([(s["name"], s["cen"], s["tel"]) for s in segs],
                         [("1", 1000, 3000), ("2", 3000, 3500), ("3", 3500, 9000)])
        self.assertEqual((js[0]["end_a"], js[0]["end_b"]), (("2", "cen"), ("3", "cen")))

    def test_ids(self):
        self.assertEqual(segment_ids(3), ["1", "2", "3"])


class FitTests(unittest.TestCase):
    """HCC1937 chr4q in miniature: copies 2 3 4 2 0 over 2 intact copies."""

    def setUp(self):
        tumor = [call((3000, "L"), (3500, "L"), 16), call((5000, "R"), (5800, "R"), 32)]
        self.js = somatic_junctions(tumor, [], WINDOW, DEFAULTS)
        self.segs = cut_segments(WINDOW, self.js)
        self.locus = Locus(self.segs, self.js)

    def test_exact_depth_gives_the_graph(self):
        unit = 15.0
        copies = {"1": 4, "2": 5, "3": 6, "4": 4, "5": 2}
        depth = {s: c * unit for s, c in copies.items()}
        best = fit(self.locus, self.js, depth, unit, unit, (Counter(), Counter()),
                   dict(DEFAULTS, starts=4))
        self.assertEqual(best["cn"], copies)
        self.assertEqual(best["mult"], {"j0": 1, "j1": 2})
        cn, mult, ref, free = remove_background(self.locus, best, (Counter(), Counter()))
        self.assertEqual(cn, {"1": 2, "2": 3, "3": 4, "4": 2, "5": 0})
        self.assertEqual(dict(free), {"1.cen": 2})

    def test_infeasible_start_is_repaired(self):
        # A bound of 3 on segment 2 forces more copies than the depth suggests; the search
        # must still find a balanced graph.
        unit = 15.0
        depth = {s: c * unit for s, c in {"1": 1, "2": 1, "3": 1, "4": 1, "5": 1}.items()}
        best = fit(self.locus, self.js, depth, unit, unit, (Counter(), Counter({"2": 3})),
                   dict(DEFAULTS, starts=2, lo=1.0, hi=1.0))
        self.assertIsNotNone(best)
        self.assertGreaterEqual(best["cn"]["2"], 3)


class GraphFileTests(unittest.TestCase):
    def graph(self):
        g = GraphFile(header={"locus": "x"})
        g.segments = [GraphSegment("1", "chr4", 0, 10, "+", 2),
                      GraphSegment("2", "chr4", 10, 20, "+", 3)]
        g.reference = [("1.R", "2.L", 2)]
        g.junctions = [("j0", "2.R", "2.R", 1, 1), ("j1", "2.L", "2.R", 1, 1)]
        g.free_ends = [("1.L", 2)]
        return g

    def test_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "graph.tsv"
            self.graph().write(path)
            again = GraphFile.read(path)
        self.assertEqual(again.segments, self.graph().segments)
        self.assertEqual(again.junctions, self.graph().junctions)
        self.assertEqual(again.header["locus"], "x")

    def test_unbalanced_is_rejected(self):
        g = self.graph()
        g.segments[1].copies = 4
        with self.assertRaises(ValueError):
            g.validate()


if __name__ == "__main__":
    unittest.main()
