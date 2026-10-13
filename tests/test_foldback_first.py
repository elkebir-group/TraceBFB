"""The fold-back-first BFB search agrees with BFB_dist on the t and total
components and with `cpp/nbfb --foldback-first` where BFB and NBFB agree."""

import subprocess
import tempfile
import unittest
from pathlib import Path

from bfbdist.distance import BFB_dist
from bfbdist.foldback_first import foldback_first_BFB_dist
from bfbdist.history import State
from bfbdist.history import Operation
from bfbdist.sequence import Sequence
from tests.test_cpp_nbfb import CASES, CPP, blocks, line


def replay(sources, ops):
    state = State(Sequence(), 0)
    for text in ops.split(";"):
        state = state.apply(Operation.from_fields(text.split(",")), sources)
    return state.sequence


class FoldbackFirstTests(unittest.TestCase):
    def test_history_spells_the_target(self):
        for strings, target in CASES:
            sources = [Sequence.from_string(line(b)) for b in blocks(strings)]
            result = foldback_first_BFB_dist(sources, Sequence.from_string(line(target)))
            with self.subTest(target=target):
                if result is None:
                    self.assertEqual(BFB_dist(sources, Sequence.from_string(line(target)))
                                     .calculate_bfb_distance(), float("inf"))
                    continue
                self.assertEqual(replay(sources, result["ops"]),
                                 Sequence.from_string(line(target)))
                self.assertEqual(len(result["types"]), result["total"])

    def test_without_foldback_joins_it_is_the_bfb_distance(self):
        # A negative span flags no join as fold-back-shaped.
        for strings, target in CASES:
            sources = [Sequence.from_string(line(b)) for b in blocks(strings)]
            target_seq = Sequence.from_string(line(target))
            result = foldback_first_BFB_dist(sources, target_seq, span=-1)
            graph = BFB_dist(sources, target_seq)
            with self.subTest(target=target):
                if result is None:
                    self.assertEqual(graph.calculate_bfb_distance(), float("inf"))
                else:
                    self.assertEqual((result["t"], result["total"]),
                                     (graph.translocation_count(),
                                      graph.calculate_bfb_distance()))

    def test_matches_cpp_nbfb(self):
        try:
            subprocess.run(["make", "-C", str(CPP)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            self.skipTest("cannot build cpp/nbfb")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "instance.txt"
            for strings, target in CASES:
                path.write_text("\n".join(map(line, strings + [target])) + "\n")
                out = subprocess.run([str(CPP / "nbfb"), "-i", str(path), "--foldback-first",
                                      "--history"], stdout=subprocess.PIPE,
                                     universal_newlines=True, check=True).stdout.split("\n")
                kv = dict(l.split(": ", 1) for l in out if ": " in l)
                total = [l for l in out if l.strip()][-1]
                sources = [Sequence.from_string(line(b)) for b in blocks(strings)]
                result = foldback_first_BFB_dist(sources, Sequence.from_string(line(target)))
                with self.subTest(target=target):
                    if result is None:
                        self.assertEqual(total, "inf")
                    else:
                        self.assertEqual((kv["foldback_t"], kv["t_cycles"], total),
                                         tuple(str(result[k])
                                               for k in ("foldback_t", "t", "total")))


if __name__ == "__main__":
    unittest.main()
