import math

from bfb_dist.bfb_graph import BFBGraph
from bfb_dist.sequence import LetterAsSegment, Segment, Sequence


def letters(text):
    return Sequence(LetterAsSegment(letter) for letter in text)


def genomic_sequence(length=3):
    return Sequence(
        Segment("chr1", index * 10, (index + 1) * 10, False)
        for index in range(length)
    )


def assert_valid_generated_graph(testcase, graph, expected_distance):
    testcase.assertIsInstance(graph, BFBGraph)
    testcase.assertEqual(graph.calculate_bfb_distance(), expected_distance)

    if expected_distance == 0 or math.isinf(expected_distance):
        testcase.assertEqual(graph.adjacency, {graph.source: {}})
        return

    testcase.assertNotIn(graph.target, graph.adjacency)
    for source, children in graph.adjacency.items():
        testcase.assertTrue(children)
        for target, operation in children.items():
            testcase.assertEqual(
                source.bfb_target(*operation),
                target,
            )

    paths = graph.get_paths()
    testcase.assertTrue(paths)
    for path in paths:
        testcase.assertEqual(path[0], graph.source)
        testcase.assertEqual(path[-1], graph.target)
        testcase.assertEqual(len(path) - 1, expected_distance)
