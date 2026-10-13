import unittest

from bfbdist.history import (
    Operation,
    State,
    start_state,
    validate_sources,
)
from bfbdist.sequence import Sequence
from tests.helpers import letters


def parse_letters(text):
    return letters(text.replace(";", ""))


class OperationTests(unittest.TestCase):
    def test_fold_back_kind_is_derived_from_r_and_flag(self):
        self.assertEqual(Operation.fold_back(0, 2, False, 2).kind, "p")
        self.assertEqual(Operation.fold_back(0, 2, True, 2).kind, "p")
        self.assertEqual(Operation.fold_back(0, 1, False, 2).kind, "i")
        self.assertEqual(Operation.fold_back(0, 1, True, 2).kind, "f")

    def test_serialization_round_trip(self):
        for operation in (
            Operation("p", 0, 2),
            Operation("i", 1, 2),
            Operation("f", 0, 1),
            Operation.translocation(3),
        ):
            with self.subTest(operation=operation):
                fields = str(operation).split(",")
                self.assertEqual(Operation.from_fields(fields), operation)

    def test_invalid_serialized_operation(self):
        for fields in (["q", "0", "1"], ["t"], ["p", "0"], ["t", "0", "1"]):
            with self.subTest(fields=fields):
                with self.assertRaises(ValueError):
                    Operation.from_fields(fields)


class StateTests(unittest.TestCase):
    def test_apply_fold_back_and_translocation(self):
        sources = (letters("AB"), letters("C"))
        state = start_state()

        state = state.apply(Operation.translocation(0), sources)
        self.assertEqual(state, State(letters("AB"), 0b01))

        state = state.apply(Operation("p", 0, 2), sources)
        self.assertEqual(state, State(letters("ABba"), 0b01))

        state = state.apply(Operation.translocation(1), sources)
        self.assertEqual(state, State(letters("ABbaC"), 2))

        state = state.apply(Operation("f", 3, 4, ), sources)
        self.assertEqual(state, State(letters("ABbacA"), 2))

        # A source may be appended again.
        state = state.apply(Operation.translocation(0), sources)
        self.assertEqual(state, State(letters("ABbacAAB"), 3))

    def test_apply_rejects_inapplicable_operations(self):
        sources = (letters("AB"),)
        state = State(letters("AB"), 0b1)

        with self.assertRaisesRegex(ValueError, "out of range"):
            state.apply(Operation.translocation(1), sources)
        with self.assertRaisesRegex(ValueError, "empty R"):
            state.apply(Operation("p", 0, 1), sources)
        with self.assertRaisesRegex(ValueError, "nonempty R"):
            state.apply(Operation("f", 0, 2), sources)
        with self.assertRaisesRegex(ValueError, "out of bounds"):
            state.apply(Operation("i", 1, 1), sources)

    def test_string_round_trip(self):
        state = State(letters("AB"), 2)

        self.assertEqual(str(state), "{2}A;B")
        self.assertEqual(State.from_string("{2}A;B", parse_letters), state)
        self.assertEqual(State.from_string("{}", parse_letters), start_state())

    def test_invalid_serialized_state(self):
        for text in ("A;B", "{0A;B", "{0,0}A;B", "{-1}A;B"):
            with self.subTest(text=text):
                with self.assertRaises(ValueError):
                    State.from_string(text, parse_letters)

    def test_translocation_count(self):
        self.assertEqual(start_state().translocations, 0)
        self.assertEqual(State(Sequence(), 2).translocations, 2)
        self.assertEqual(str(State(Sequence(), 2)), "{2}")

    def test_validate_sources_rejects_empty_source(self):
        with self.assertRaisesRegex(ValueError, "source 1 is empty"):
            validate_sources([letters("A"), letters("")])


if __name__ == "__main__":
    unittest.main()
