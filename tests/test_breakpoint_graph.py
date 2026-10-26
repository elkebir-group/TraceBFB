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

    def test_lwalk_key_leaves_out_reference_adjacencies(self):
        # 60 segments of 30 copies: the reference adjacencies alone would need 295 key bits,
        # but they follow from the one fold-back, so the walk search fits
        try:
            subprocess.run(["make", "-C", str(CPP)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            self.skipTest("cannot build cpp/lwalk")
        n = 60
        rows = ["H\tversion=1\tlocus=chain"]
        rows += [f"S\t{i}\tchr1\t{100 * i}\t{100 * i + 100}\t+\t30" for i in range(1, n + 1)]
        rows += [f"R\t{i}.R\t{i + 1}.L\t30" for i in range(1, n)]
        rows += [f"J\tj0\t{n}.R\t{n}.R\t15\t1", "F\t1.L\t30"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "graph.tsv"
            path.write_text("\n".join(rows) + "\n")
            out = subprocess.run([str(CPP / "lwalk"), str(path), "--max-foldbacks"],
                                 capture_output=True, text=True, check=True).stdout
        self.assertIn(f"best walk: {2 * n} of {30 * n} traversals, 1 fold-back traversals", out)
        self.assertIn("walks of that value 1", out)

    def test_lwalk_too_large_exits_5(self):
        try:
            subprocess.run(["make", "-C", str(CPP)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            self.skipTest("cannot build cpp/lwalk")
        # one segment of 20001 copies in tandem: over the 20000-traversal recursion limit
        rows = ["H\tversion=1\tlocus=tandem", "S\t1\tchr1\t0\t100\t+\t20001",
                "J\tj0\t1.R\t1.L\t20000\t0", "F\t1.L\t1", "F\t1.R\t1"]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "graph.tsv"
            path.write_text("\n".join(rows) + "\n")
            run = subprocess.run([str(CPP / "lwalk"), str(path)], capture_output=True, text=True)
        self.assertEqual(run.returncode, 5)
        self.assertIn("Too large", run.stderr)


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
