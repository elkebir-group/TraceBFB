"""The C++ general-source search (`cpp/nbfb --simple`) agrees with the C++
string search under the same definition and with a brute-force search that
applies the foldable-prefix definition literally. With `--blocks`, the
search over all blocks of the source lines agrees with the general search
over the explicit list of blocks. With `--all-optimal`, both list exactly
the optimal histories of a brute-force enumeration (a history is its
strings and cycle types)."""

import heapq
import random
from collections import defaultdict
import subprocess
import tempfile
import unittest
from pathlib import Path

CPP = Path(__file__).resolve().parents[1] / "cpp"

# Letters A-C lie on chr1 and D on chr2; upper case is forward, lower case
# reverse. A join is fold-back-shaped when both segments lie on the same
# chromosome in opposite orientations (all segments are within 50 kb).
CHROMOSOME = {"A": "chr1", "B": "chr1", "C": "chr1", "D": "chr2"}


def rc(s):
    return s[::-1].swapcase()


def token(letter):
    k = 10 * (ord(letter.upper()) - 64)
    chromosome = CHROMOSOME[letter.upper()]
    return (f"{chromosome}:{k}-{k + 9}[R]" if letter.isupper()
            else f"{chromosome}:{k + 9}-{k}[C]")


def line(text):
    return ";".join(token(c) for c in text)


def foldable(s, t):
    """s = x y with x a prefix of t and y or rc(y) a substring of t."""
    return any(t.startswith(s[:a]) and (s[a:] in t or rc(s[a:]) in t)
               for a in range(len(s) + 1))


def foldback(u, v):
    return CHROMOSOME[u.upper()] == CHROMOSOME[v.upper()] and u.isupper() != v.isupper()


def brute(sources, t, foldback_first):
    """Minimum (fold-back t, t, total) over nested histories, or None."""
    best = {"": (0, 0, 0)}
    heap = [((0, 0, 0), "")]
    while heap:
        cost, s = heapq.heappop(heap)
        if best[s] != cost:
            continue
        if s == t:
            return cost
        children = []
        for a in range(len(s)):
            for b in range(a + 1, len(s) + 1):
                children.append((s + rc(s[a:b]), (0, 0, 1)))
                if b < len(s):
                    children.append((s[:b] + rc(s[b:]) + rc(s[a:b]), (0, 0, 1)))
        for x in sources:
            flag = int(foldback_first and len(s) > 0 and foldback(s[-1], x[0]))
            children.append((s + x, (flag, 1, 1)))
        for w, step in children:
            if len(w) > len(t) or not foldable(w, t):
                continue
            new = tuple(c + d for c, d in zip(cost, step))
            if w not in best or new < best[w]:
                best[w] = new
                heapq.heappush(heap, (new, w))
    return None


def children(s, sources, foldback_first):
    """(string, cycle type, cost step) of every cycle on s."""
    for a in range(len(s)):
        for b in range(a + 1, len(s) + 1):
            q = rc(s[a:b])
            yield s + q, "p" if b == len(s) else "i", (0, 0, 1)
            if b < len(s):
                yield s[:b] + rc(s[b:]) + q, "f", (0, 0, 1)
    for x in sources:
        flag = int(foldback_first and len(s) > 0 and foldback(s[-1], x[0]))
        yield s + x, "t", (flag, 1, 1)


def all_optimal(sources, t, foldback_first):
    """Every optimal nested history, as (strings, cycle types)."""
    best = {"": (0, 0, 0)}
    heap = [((0, 0, 0), "")]
    done = set()
    while heap:
        cost, s = heapq.heappop(heap)
        if s in done or best[s] != cost:
            continue
        done.add(s)
        for w, _, step in children(s, sources, foldback_first):
            if len(w) > len(t) or not foldable(w, t):
                continue
            new = tuple(c + d for c, d in zip(cost, step))
            if w not in best or new < best[w]:
                best[w] = new
                heapq.heappush(heap, (new, w))
    if t not in best:
        return set()
    into = defaultdict(set)
    for s in done:
        for w, kind, step in children(s, sources, foldback_first):
            if w in best and tuple(c + d for c, d in zip(best[s], step)) == best[w]:
                into[w].add((s, kind))
    histories = set()

    def walk(v, strings, kinds):
        if v == "":
            histories.add((tuple(reversed(strings)), tuple(reversed(kinds))))
            return
        for u, kind in into[v]:
            walk(u, strings + [u], kinds + [kind])

    walk(t, [t], [])
    return histories


def replay(ops, sources):
    """(strings, cycle types) of a history printed by nbfb."""
    s, strings, kinds = "", [""], []
    for op in ops.split(";"):
        kind, *args = op.split(",")
        if kind == "t":
            s += sources[int(args[0])]
        else:
            x, y = int(args[0]), int(args[1])
            p, q, r = s[:x], s[x:y], s[y:]
            s = p + q + r + rc(q) if kind in "pi" else p + q + rc(r) + rc(q)
        strings.append(s)
        kinds.append(kind)
    return tuple(strings), tuple(kinds)


def blocks(strings):
    """Every block of every string, in both orientations, once each."""
    out = []
    for s in strings:
        for i in range(len(s)):
            for j in range(i + 1, len(s) + 1):
                for b in (s[i:j], rc(s[i:j])):
                    if b not in out:
                        out.append(b)
    return out


def random_instance(rng):
    alphabet = "ABCDabcd"
    sources = ["".join(rng.choice(alphabet) for _ in range(rng.randint(1, 3)))
               for _ in range(rng.randint(1, 2))]
    t = ""
    while len(t) < rng.randint(2, 7):
        if not t or rng.random() < 0.4:
            t += rng.choice(sources)
        else:
            b = rng.randint(1, len(t))
            a = rng.randrange(b)
            if b == len(t) or rng.random() < 0.5:
                t += rc(t[a:b])
            else:
                t = t[:b] + rc(t[b:]) + rc(t[a:b])
    if rng.random() < 0.3:
        sources.append("".join(rng.choice(alphabet) for _ in range(2)))
    return list(dict.fromkeys(sources)), t


class CppSimpleSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            subprocess.run(["make", "-C", str(CPP)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except (OSError, subprocess.CalledProcessError):
            raise unittest.SkipTest("cannot build cpp/nbfb")

    def run_nbfb(self, path, *flags):
        result = subprocess.run(
            [str(CPP / "nbfb"), "-i", str(path), "--history", "--max-mem-gb", "1", *flags],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, check=True)
        lines = [l for l in result.stdout.split("\n") if l.strip()]
        fields = dict(l.split(": ", 1) for l in lines[:-1] if ": " in l)
        if lines[-1] == "inf":
            return None
        return (int(fields.get("foldback_t", 0)), int(fields["t_cycles"]), int(lines[-1]))

    def test_simple_search_matches_brute_force_and_string_search(self):
        rng = random.Random(7)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "instance.txt"
            for _ in range(80):
                sources, t = random_instance(rng)
                path.write_text("\n".join(map(line, sources + [t])) + "\n")
                for foldback_first in (False, True):
                    flags = ["--foldback-first"] if foldback_first else []
                    with self.subTest(sources=sources, t=t, flags=flags):
                        expected = brute(sources, t, foldback_first)
                        if expected is not None and not foldback_first:
                            expected = (0,) + expected[1:]
                        self.assertEqual(self.run_nbfb(path, "--simple", *flags), expected)
                        self.assertEqual(
                            self.run_nbfb(path, "--simple", "--string-search", *flags),
                            expected)

    def test_block_search_matches_brute_force_over_all_blocks(self):
        rng = random.Random(11)
        with tempfile.TemporaryDirectory() as directory:
            given = Path(directory) / "given.txt"
            closed = Path(directory) / "closed.txt"
            for _ in range(60):
                strings, t = random_instance(rng)
                given.write_text("\n".join(map(line, strings + [t])) + "\n")
                closed.write_text("\n".join(map(line, blocks(strings) + [t])) + "\n")
                for foldback_first in (False, True):
                    flags = ["--foldback-first"] if foldback_first else []
                    with self.subTest(strings=strings, t=t, flags=flags):
                        expected = brute(blocks(strings), t, foldback_first)
                        if expected is not None and not foldback_first:
                            expected = (0,) + expected[1:]
                        self.assertEqual(
                            self.run_nbfb(given, "--simple", "--blocks", *flags), expected)
                        self.assertEqual(self.run_nbfb(closed, "--simple", *flags), expected)

    def run_all_optimal(self, path, sources, *flags):
        result = subprocess.run(
            [str(CPP / "nbfb"), "-i", str(path), "--simple", "--all-optimal",
             "--max-mem-gb", "1", *flags],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, check=True)
        lines = [l for l in result.stdout.split("\n") if l.strip()]
        count = int(next(l for l in lines if l.startswith("optimal_histories:")).split(": ")[1])
        listed = [replay(l.split("\t")[1], sources) for l in lines if l.startswith("optimal: ")]
        return count, listed

    def test_all_optimal_histories_match_brute_force(self):
        rng = random.Random(5)
        with tempfile.TemporaryDirectory() as directory:
            given = Path(directory) / "given.txt"
            closed = Path(directory) / "closed.txt"
            for _ in range(40):
                strings, t = random_instance(rng)
                every = blocks(strings)
                given.write_text("\n".join(map(line, strings + [t])) + "\n")
                closed.write_text("\n".join(map(line, every + [t])) + "\n")
                for foldback_first in (False, True):
                    flags = ["--foldback-first"] if foldback_first else []
                    with self.subTest(strings=strings, t=t, flags=flags):
                        expected = all_optimal(strings, t, foldback_first)
                        count, listed = self.run_all_optimal(given, strings, *flags)
                        self.assertEqual(count, len(expected))
                        self.assertEqual(len(listed), len(set(listed)))
                        self.assertEqual(set(listed), expected)
                        # Block search, every block listed so that each appended
                        # block is named.
                        expected = all_optimal(every, t, foldback_first)
                        count, listed = self.run_all_optimal(closed, every, "--blocks", *flags)
                        self.assertEqual(count, len(expected))
                        self.assertEqual(set(listed), expected)


if __name__ == "__main__":
    unittest.main()
