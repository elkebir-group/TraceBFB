from contextlib import redirect_stdout
from io import StringIO
import re
import unittest

from tracebfb.history import start_state
from tracebfb.sequence import Sequence
from tracebfb.simulation import (
    _set_seed,
    simulate_history,
    simulate_reference_sequence,
    simulate_sources,
)
from tests.helpers import genomic_sequence, letters


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
            (2, 10, 20000, None, "chrM"),
        )

        for arguments in invalid_arguments:
            with self.subTest(arguments=arguments):
                with self.assertRaises(ValueError):
                    simulate_reference_sequence(*arguments)

    def test_sources_lie_on_distinct_chromosomes(self):
        _set_seed(5)
        sources = simulate_sources(3, 2, 10, 20)

        self.assertEqual(len(sources), 3)
        chromosomes = {source[0].chr for source in sources}
        self.assertEqual(len(chromosomes), 3)
        for source in sources:
            self.assertEqual(len(source), 2)
            self.assertEqual({segment.chr for segment in source.segments}, {source[0].chr})

    def test_invalid_source_count(self):
        with self.assertRaises(ValueError):
            simulate_sources(0, 2)


class HistorySimulationTests(unittest.TestCase):
    def assert_valid_history(self, steps, sources, num_cycles):
        self.assertEqual(len(steps), num_cycles + 1)
        self.assertEqual(steps[0][0], start_state())
        self.assertEqual(steps[0][1].kind, "t")
        self.assertIsNone(steps[-1][1])
        appended = []
        for (state, operation), (child, _) in zip(steps, steps[1:]):
            self.assertEqual(state.apply(operation, sources), child)
            if operation.kind == "t":
                appended.append(operation.source)
        self.assertEqual(sorted(appended), list(range(len(sources))))

    def test_history_is_reproducible_and_silent(self):
        sources = [genomic_sequence(3, "chr1"), genomic_sequence(2, "chr2")]
        stdout = StringIO()

        _set_seed(29)
        with redirect_stdout(stdout):
            first_steps, first_color = simulate_history(sources, 5, False, True)
        _set_seed(29)
        second_steps, second_color = simulate_history(sources, 5, False, True)

        self.assertEqual(stdout.getvalue(), "")
        self.assertEqual(first_steps, second_steps)
        self.assertEqual(first_color, second_color)
        self.assert_valid_history(first_steps, sources, 5)
        self.assertEqual(
            [ANSI.sub("", line) for line in first_color.splitlines()],
            [str(state.sequence) for state, _ in first_steps[1:]],
        )

    def test_short_deletions_bound_every_forward_overhang(self):
        sources = [letters("ABCD"), letters("EF")]
        for seed in range(30):
            _set_seed(seed)
            steps, color_output = simulate_history(sources, 7, True)
            self.assertIsNone(color_output)
            self.assert_valid_history(steps, sources, 7)
            previous_gain = None
            for state, operation in steps[:-1]:
                if operation.kind == "t":
                    gain = len(sources[operation.source])
                else:
                    gain = operation.y - operation.x
                    if operation.kind == "f":
                        overhang = len(state.sequence) - operation.y
                        self.assertLess(overhang, previous_gain)
                previous_gain = gain

    def test_single_source_single_cycle(self):
        source = genomic_sequence(2)

        steps, color_output = simulate_history([source], 1, False, True)

        self.assertEqual(len(steps), 2)
        self.assertEqual(steps[-1][0].sequence, source)
        self.assertEqual(color_output, f"{source}\n")

    def test_invalid_history_arguments(self):
        with self.assertRaisesRegex(ValueError, "at least one source"):
            simulate_history([], 1)
        with self.assertRaisesRegex(ValueError, "cannot be empty"):
            simulate_history([Sequence()], 1)
        with self.assertRaisesRegex(ValueError, "cannot be smaller"):
            simulate_history([genomic_sequence(1)] * 2, 1)


if __name__ == "__main__":
    unittest.main()
