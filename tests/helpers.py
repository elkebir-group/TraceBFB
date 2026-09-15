import math

from bfb_dist.bfb_graph import BFBGraph
from bfb_dist.sequence import LetterAsSegment, Segment, Sequence


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

    testcase.assertNotIn(graph.end, graph.adjacency)
    for node, children in graph.adjacency.items():
        testcase.assertTrue(children)
        for child, operation in children.items():
            testcase.assertEqual(node.apply(operation, graph.sources), child)

    paths = graph.get_paths()
    testcase.assertTrue(paths)
    for path in paths:
        testcase.assertEqual(path[0], graph.start)
        testcase.assertEqual(path[-1], graph.end)
        testcase.assertEqual(len(path) - 1, expected_distance)
        kinds = [
            graph.adjacency[node][child].kind
            for node, child in zip(path, path[1:])
        ]
        testcase.assertEqual(kinds[0], "t")
        testcase.assertEqual(kinds.count("t"), len(graph.sources))
