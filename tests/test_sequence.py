import unittest

from bfb_dist.sequence import Segment, Sequence
from tests.helpers import genomic_sequence, letters


class SegmentTests(unittest.TestCase):
    def test_string_and_reverse_complement(self):
        segment = Segment("chr2", 10, 25, False)

        self.assertEqual(str(segment), "chr2:10-25[R]")
        self.assertEqual(
            segment.reverse_complement(),
            Segment("chr2", 25, 10, True),
        )
        self.assertEqual(
            segment.reverse_complement().reverse_complement(),
            segment,
        )


class SequenceTests(unittest.TestCase):
    def test_hashing_slicing_and_addition(self):
        sequence = genomic_sequence(3)

        self.assertEqual(sequence[:2], Sequence(sequence.segments[:2]))
        self.assertEqual(sequence[0], sequence.segments[0])
        self.assertEqual(sequence[:1] + sequence[1:], sequence)
        self.assertEqual(sequence[:2] + sequence[2], sequence)
        self.assertEqual({sequence: "value"}[sequence], "value")

    def test_decompose_supports_empty_p_and_r(self):
        sequence = letters("ABC")

        p, q, r = sequence.decompose(0, 2)
        self.assertEqual((str(p), str(q), str(r)), ("", "A;B", "C"))

        p, q, r = sequence.decompose(1, len(sequence))
        self.assertEqual((str(p), str(q), str(r)), ("A", "B;C", ""))

    def test_bfb_target_formulas(self):
        sequence = letters("ABC")

        self.assertEqual(
            sequence.bfb_target(1, 2, False),
            letters("ABCb"),
        )
        self.assertEqual(
            sequence.bfb_target(1, 2, True),
            letters("ABcb"),
        )

    def test_empty_r_makes_deletion_direction_equivalent(self):
        sequence = letters("ABC")

        self.assertEqual(
            sequence.bfb_target(1, 3, False),
            sequence.bfb_target(1, 3, True),
        )

    def test_previous_bfb_recovers_source_and_operation(self):
        source = letters("ABC")
        target = source.bfb_target(1, 2, True)

        self.assertIn((source, 1, 2, True), target.previous_bfb())

    def test_genomic_string_round_trip(self):
        sequence = genomic_sequence(3)

        self.assertEqual(Sequence.from_string(str(sequence)), sequence)

    def test_invalid_serialized_segment_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Invalid segment"):
            Sequence.from_string("chr1:10-20[X]")


if __name__ == "__main__":
    unittest.main()
