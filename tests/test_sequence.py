import unittest

from bfbdist.sequence import Segment, Sequence
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

    def test_successors_produces_forward_and_non_forward_children(self):
        source = letters("ABC")

        successors = source.successors()

        self.assertIn(
            (source.bfb_target(1, 2, False), 1, 2, False),
            successors,
        )
        self.assertIn(
            (source.bfb_target(1, 2, True), 1, 2, True),
            successors,
        )

    def test_successors_omits_forward_deletion_when_r_is_empty(self):
        source = letters("ABC")

        successors = source.successors()

        self.assertNotIn(
            (source.bfb_target(0, 3, True), 0, 3, True),
            successors,
        )
        self.assertIn(
            (source.bfb_target(0, 3, False), 0, 3, False),
            successors,
        )

    def test_successors_and_previous_bfb_are_inverse(self):
        source = letters("ABC")

        for child, x, y, is_forward_deletion in source.successors():
            self.assertIn(
                (source, x, y, is_forward_deletion),
                child.previous_bfb(),
            )

    def test_suffix_helpers(self):
        self.assertTrue(letters("BC").is_suffix_of(letters("ABC")))
        self.assertTrue(letters("").is_suffix_of(letters("ABC")))
        self.assertFalse(letters("AB").is_suffix_of(letters("ABC")))
        self.assertFalse(letters("ABCD").is_suffix_of(letters("ABC")))
        self.assertEqual(letters("ABC").strip_suffix(letters("BC")), letters("A"))
        self.assertEqual(letters("ABC").strip_suffix(letters("ABC")), letters(""))
        self.assertIsNone(letters("ABC").strip_suffix(letters("AB")))

    def test_lcp_length(self):
        self.assertEqual(
            letters("ABC").lcp_length(letters("ABD")),
            2,
        )
        self.assertEqual(
            letters("ABC").lcp_length(letters("ABC")),
            3,
        )
        self.assertEqual(
            letters("ABC").lcp_length(letters("")),
            0,
        )

    def test_is_infix_of(self):
        self.assertTrue(letters("BC").is_infix_of(letters("ABCD")))
        self.assertTrue(letters("").is_infix_of(letters("ABCD")))
        self.assertFalse(letters("CB").is_infix_of(letters("ABCD")))
        self.assertFalse(letters("BC").is_infix_of(letters("AB")))

    def test_is_infix_of_boundary_matches(self):
        # Match at the very start and the very end of the haystack.
        self.assertTrue(letters("AB").is_infix_of(letters("ABCD")))
        self.assertTrue(letters("CD").is_infix_of(letters("ABCD")))
        self.assertTrue(letters("ABCD").is_infix_of(letters("ABCD")))
        self.assertFalse(letters("ABCDE").is_infix_of(letters("ABCD")))

    def test_is_infix_of_overlapping_repeats(self):
        # Classic Z-algorithm stress cases: repeated/overlapping patterns.
        self.assertTrue(letters("AA").is_infix_of(letters("AAAA")))
        self.assertTrue(letters("ABAB").is_infix_of(letters("ABABAB")))
        self.assertFalse(letters("AAB").is_infix_of(letters("AAAA")))
        self.assertTrue(letters("AAB").is_infix_of(letters("AAAAB")))

    def test_genomic_string_round_trip(self):
        sequence = genomic_sequence(3)

        self.assertEqual(Sequence.from_string(str(sequence)), sequence)

    def test_invalid_serialized_segment_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Invalid segment"):
            Sequence.from_string("chr1:10-20[X]")


if __name__ == "__main__":
    unittest.main()
