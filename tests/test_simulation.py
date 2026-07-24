from contextlib import redirect_stdout
from io import StringIO
import re
import unittest

from bfb_dist.sequence import Sequence
from bfb_dist.simulation import (
    _set_seed,
    simulate_bfb,
    simulate_reference_sequence,
)
from tests.helpers import genomic_sequence


ANSI = re.compile(r"\033\[[0-9;]*m")


class ReferenceSimulationTests(unittest.TestCase):
    def test_reference_is_reproducible_and_contiguous(self):
        _set_seed(17)
        first = simulate_reference_sequence(4, 10, 10, 100)
        _set_seed(17)
        second = simulate_reference_sequence(4, 10, 10, 100)

        self.assertEqual(first, second)
        self.assertEqual(len(first), 4)
        self.assertEqual(first[0].start, 100)
        for left, right in zip(first.segments, first.segments[1:]):
            self.assertEqual(left.end, right.start)

    def test_invalid_reference_arguments(self):
        invalid_arguments = (
            (0, 10, 20, None),
            (2, 0, 20, None),
            (2, 20, 10, None),
            (2, 10, 20, -1),
        )

        for arguments in invalid_arguments:
            with self.subTest(arguments=arguments):
                with self.assertRaises(ValueError):
                    simulate_reference_sequence(*arguments)


class BFBSimulationTests(unittest.TestCase):
    def test_simulation_is_reproducible_and_silent(self):
        reference = genomic_sequence(3)
        stdout = StringIO()

        _set_seed(29)
        with redirect_stdout(stdout):
            first_ops, first_color = simulate_bfb(reference, 3, True)
        _set_seed(29)
        second_ops, second_color = simulate_bfb(reference, 3, True)

        self.assertEqual(stdout.getvalue(), "")
        self.assertEqual(first_ops, second_ops)
        self.assertEqual(first_color, second_color)
        self.assertEqual(len(first_ops), 4)
        self.assertIn("\033[91m", first_color)
        self.assertEqual(
            [
                ANSI.sub("", line)
                for line in first_color.splitlines()
            ],
            list(map(str, first_ops)),
        )

        for source, target in zip(first_ops, first_ops[1:]):
            self.assertTrue(
                any(
                    predecessor == source
                    for predecessor, _, _, _ in target.previous_bfb()
                )
            )

    def test_uncolored_simulation_returns_no_color_string(self):
        operations, color_output = simulate_bfb(
            genomic_sequence(2),
            0,
            False,
        )

        self.assertEqual(len(operations), 1)
        self.assertIsNone(color_output)

    def test_colored_zero_cycle_representation(self):
        reference = genomic_sequence(2)

        operations, color_output = simulate_bfb(reference, 0, True)

        self.assertEqual(operations, [reference])
        self.assertEqual(color_output, f"{reference}\n")

    def test_invalid_cycle_arguments(self):
        with self.assertRaisesRegex(ValueError, "cannot be negative"):
            simulate_bfb(genomic_sequence(1), -1)
        with self.assertRaisesRegex(ValueError, "cannot be empty"):
            simulate_bfb(Sequence(), 1)


if __name__ == "__main__":
    unittest.main()
