"""The C++ block search (`cpp/nbfb`) agrees with the C++ string search and
with NBFB_dist when the sources are all blocks of the source lines."""

import subprocess
import tempfile
import unittest
from pathlib import Path

from bfbdist.nbfb import NBFB_dist
from bfbdist.sequence import Sequence

CPP = Path(__file__).resolve().parents[1] / "cpp"

# (source strings, target); upper case is forward, lower case reverse.
CASES = [
    (["ABCD"], "ABCDdcb"),
    (["ABCD", "X"], "ABCDXcb"),
    (["ABCD", "x"], "ABCXdcb"),
    (["ABCD", "y"], "ABCDcYBC"),
    (["AAba"], "AAAAAABaa"),
    (["AAB"], "AABbAbAAB"),
    (["ZABCDEFGP"], "ZABCFEFGPgfefcbaB"),
    (["AB"], "C"),
]


def token(letter):
    k = 10 * (ord(letter.upper()) - 64)
    return f"chr1:{k}-{k + 9}[R]" if letter.isupper() else f"chr1:{k + 9}-{k}[C]"


def line(text):
    return ";".join(token(c) for c in text)


def blocks(strings):
    """Every block of every string, in both orientations, once each."""
    out = []
    for s in strings:
        for i in range(len(s)):
            for j in range(i + 1, len(s) + 1):
                for b in (s[i:j], s[i:j][::-1].swapcase()):
                    if b not in out:
                        out.append(b)
    return out


class CppBlockSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            subprocess.run(["make", "-C", str(CPP)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            raise unittest.SkipTest("cannot build cpp/nbfb")

    def run_nbfb(self, path, *flags):
        result = subprocess.run([str(CPP / "nbfb"), "-i", str(path), "--history", *flags],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, check=True)
        lines = result.stdout.split("\n")
        t_cycles = next(l for l in lines if l.startswith("t_cycles:")).split(": ")[1]
        return t_cycles, [l for l in lines if l.strip()][-1]

    def test_block_search_matches_string_search(self):
        with tempfile.TemporaryDirectory() as directory:
            given, closed = Path(directory) / "given.txt", Path(directory) / "closed.txt"
            for strings, target in CASES:
                given.write_text("\n".join(map(line, strings + [target])) + "\n")
                closed.write_text("\n".join(map(line, blocks(strings) + [target])) + "\n")
                for flags in ([], ["--foldback-first"]):
                    with self.subTest(strings=strings, target=target, flags=flags):
                        self.assertEqual(self.run_nbfb(given, *flags),
                                         self.run_nbfb(closed, "--string-search", *flags))
                with self.subTest(strings=strings, target=target, python=True):
                    graph = NBFB_dist([Sequence.from_string(line(b)) for b in blocks(strings)],
                                      Sequence.from_string(line(target)))
                    distance = graph.calculate_bfb_distance()
                    expected = ("inf", "inf") if distance == float("inf") else (
                        str(graph.translocation_count()), str(distance))
                    self.assertEqual(self.run_nbfb(given), expected)


if __name__ == "__main__":
    unittest.main()
