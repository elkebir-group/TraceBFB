"""Breakpoint graphs of assembled targets and their walks by `cpp/lwalk`."""

import subprocess
import tempfile
import unittest
from pathlib import Path

from bfbdist.breakpoint_graph import BreakpointGraph
from bfbdist.sequence import Sequence
from tests.helpers import CPP, line

# HG008-T chr19 (Z A B C D E F G P, D deleted) and a single fold-back.
CASES = [
    ("ZABCDEFGP", "ZABCFEFGPgfefcbaB", 12),
    ("ABC", "ABCcb", 1),
]


class BreakpointGraphTests(unittest.TestCase):
    def test_reference_multiplicities_follow_from_balance(self):
        graph = BreakpointGraph.from_target([Sequence.from_string(line("ABC"))],
                                            Sequence.from_string(line("ABCcb")))
        self.assertEqual(sorted(graph.copy_numbers.values()), [1, 2, 2])
        self.assertEqual(sum(graph.junctions.values()), 1)
        self.assertEqual(sum(graph.reference.values()), 3)

    def test_lwalk_enumerates_the_walks(self):
        try:
            subprocess.run(["make", "-C", str(CPP)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            self.skipTest("cannot build cpp/lwalk")
        with tempfile.TemporaryDirectory() as directory:
            for chain, target, count in CASES:
                target_seq = Sequence.from_string(line(target))
                graph = BreakpointGraph.from_target([Sequence.from_string(line(chain))],
                                                    target_seq)
                path, prefix = Path(directory) / "graph.txt", Path(directory) / "walks"
                graph.write_lwalk(path)
                subprocess.run([str(CPP / "lwalk"), str(path), "--enumerate", "-o", str(prefix)],
                               stdout=subprocess.DEVNULL, check=True)
                walks = Path(f"{prefix}.targets.txt").read_text().split()
                with self.subTest(target=target):
                    self.assertEqual(len(walks), count)
                    self.assertIn(str(target_seq), walks)
                    lengths = {len(w.split(";")) for w in walks}
                    self.assertEqual(lengths, {len(target)})
                    if target == "ABCcb":
                        # default names s0, s1, ... have no lower case: reverse is -name
                        self.assertEqual(Path(f"{prefix}.walks.txt").read_text().split(),
                                         ["s0", "s1", "s2", "-s2", "-s1"])


class FoldbackTests(unittest.TestCase):
    def test_reverse_read_segment_uses_genomic_sides(self):
        from bfbdist.breakpoint_graph import _foldback
        fwd = Sequence.from_string("chr1:100-200[R]")[0]
        rev = Sequence.from_string("chr1:400-300[R]")[0]  # read right to left: start 400 > end 300
        # end of fwd (genomic R, 200) and start of rev (genomic R, 400): a fold-back
        self.assertTrue(_foldback(frozenset({(fwd, "t"), (rev, "h")}), 50_000))
        # end of fwd (R) and end of rev (genomic L, 300): not a fold-back
        self.assertFalse(_foldback(frozenset({(fwd, "t"), (rev, "t")}), 50_000))


if __name__ == "__main__":
    unittest.main()
