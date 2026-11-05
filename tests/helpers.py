from pathlib import Path
import math

from tracebfb.bfb_graph import BFBGraph
from tracebfb.sequence import LetterAsSegment, Segment, Sequence


def letters(text):
    return Sequence(LetterAsSegment(letter) for letter in text)


def genomic_sequence(length=3, chromosome="chr1"):
    return Sequence(
        Segment(chromosome, index * 10, (index + 1) * 10, False)
        for index in range(length)
    )


def assert_valid_generated_graph(testcase, graph, expected_distance):
    testcase.assertIsInstance(graph, BFBGraph)
    testcase.assertEqual(graph.calculate_bfb_distance(), expected_distance)

    if expected_distance == 0 or math.isinf(expected_distance):
        testcase.assertEqual(graph.adjacency, {graph.start: {}})
        return

    for node in graph.adjacency:
        testcase.assertFalse(graph.is_end(node))
    for node, children in graph.adjacency.items():
        testcase.assertTrue(children)
        for child, operation in children.items():
            testcase.assertEqual(node.apply(operation, graph.sources), child)

    paths = graph.get_paths()
    testcase.assertTrue(paths)
    for path in paths:
        testcase.assertEqual(path[0], graph.start)
        testcase.assertTrue(graph.is_end(path[-1]))
        testcase.assertEqual(len(path) - 1, expected_distance)
        kinds = [
            graph.adjacency[node][child].kind
            for node, child in zip(path, path[1:])
        ]
        testcase.assertEqual(kinds[0], "t")
        testcase.assertEqual(kinds.count("t"), path[-1].translocations)


# Letter instances for the C++ tools: upper case is a forward segment on chr1,
# lower case its reverse complement.
CPP = Path(__file__).resolve().parents[1] / "cpp"

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
