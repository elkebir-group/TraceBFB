import unittest

from bfb_dist.utils import format_sequence
from tests.helpers import letters


class FormattingTests(unittest.TestCase):
    def test_plain_sequence(self):
        parts = (letters("A"), letters("B"), letters("C"))

        self.assertEqual(format_sequence(parts), "A;B;C")

    def test_colored_sequence(self):
        parts = (letters("A"), letters("B"), letters("C"))

        formatted = format_sequence(parts, colored=True)

        self.assertIn("\033[91mA\033[00m", formatted)
        self.assertIn("\033[92mB\033[00m", formatted)
        self.assertIn("\033[95mC\033[00m", formatted)
        self.assertEqual(formatted.count(";"), 2)

    def test_empty_p_and_r_do_not_add_separators(self):
        parts = (letters(""), letters("B"), letters(""))

        self.assertEqual(format_sequence(parts), "B")
        self.assertEqual(
            format_sequence(parts, colored=True),
            "\033[92mB\033[00m",
        )


if __name__ == "__main__":
    unittest.main()
