"""States and operations of BFB histories over a set of sources.

A history starts from the empty string and applies p-, i-, f-, and t-BFB
cycles. A t-BFB cycle appends one source from the input set, and every
source is appended exactly once. A search state therefore pairs the current
sequence with the set of sources appended so far, stored as a bitmask over
source indices so that equal sources at different indices stay distinct.
"""

from dataclasses import dataclass

from bfb_dist.sequence import Sequence


KINDS = ("p", "i", "f", "t")


@dataclass(frozen=True)
class Operation:
    """One BFB cycle.

    For kinds ``p``, ``i``, and ``f``, the decomposition of the sequence the
    cycle acts on is ``P = v[:x]``, ``Q = v[x:y]``, and ``R = v[y:]``. Kind
    ``p`` requires ``R`` to be empty, kinds ``i`` and ``f`` require ``R`` to
    be nonempty. For kind ``t``, ``source`` is the index of the appended
    source in the source list.
    """

    kind: str
    x: int = -1
    y: int = -1
    source: int = -1

    @classmethod
    def fold_back(cls, x, y, is_forward_deletion, length):
        """Build a p-, i-, or f-BFB operation from the legacy triple."""
        if y >= length:
            return cls("p", x, y)
        if is_forward_deletion:
            return cls("f", x, y)
        return cls("i", x, y)

    @classmethod
    def translocation(cls, source):
        """Build a t-BFB operation appending source ``source``."""
        return cls("t", source=source)

    @property
    def is_forward_deletion(self):
        """Return whether the operation deletes the forward overhang."""
        return self.kind == "f"

    def __str__(self):
        """Serialize as ``kind,x,y`` or ``t,source``."""
        if self.kind == "t":
            return f"t,{self.source}"
        return f"{self.kind},{self.x},{self.y}"

    @classmethod
    def from_fields(cls, fields):
        """Parse the fields of a serialized operation.

        Raises:
            ValueError: If the fields do not describe an operation.
        """
        if len(fields) == 2 and fields[0] == "t":
            return cls.translocation(int(fields[1]))
        if len(fields) == 3 and fields[0] in ("p", "i", "f"):
            return cls(fields[0], int(fields[1]), int(fields[2]))
        raise ValueError(f"invalid operation {','.join(fields)!r}")


def full_mask(sources):
    """Return the bitmask with every source index set."""
    return (1 << len(sources)) - 1


def mask_indices(mask):
    """Return the sorted source indices contained in ``mask``."""
    indices = []
    index = 0
    while mask:
        if mask & 1:
            indices.append(index)
        mask >>= 1
        index += 1
    return indices


@dataclass(frozen=True)
class State:
    """A sequence together with the bitmask of sources appended so far."""

    sequence: Sequence
    used: int = 0

    def __str__(self):
        """Serialize as ``{i,j,...}`` followed by the sequence."""
        indices = ",".join(map(str, mask_indices(self.used)))
        return "{" + indices + "}" + str(self.sequence)

    __repr__ = __str__

    def __len__(self):
        """Return the number of segments of the sequence."""
        return len(self.sequence)

    def has_used(self, index):
        """Return whether source ``index`` has been appended."""
        return bool(self.used >> index & 1)

    def apply(self, operation, sources):
        """Return the state produced by ``operation``.

        Raises:
            ValueError: If the operation is not applicable to this state.
        """
        sequence = self.sequence
        if operation.kind == "t":
            index = operation.source
            if not 0 <= index < len(sources):
                raise ValueError(
                    f"source index {index} out of range"
                )
            if self.has_used(index):
                raise ValueError(
                    f"source {index} has already been appended"
                )
            return State(sequence + sources[index], self.used | 1 << index)

        x, y = operation.x, operation.y
        length = len(sequence)
        if not 0 <= x < y <= length:
            raise ValueError(
                f"decomposition ({x}, {y}) is out of bounds for a "
                f"sequence of length {length}"
            )
        if operation.kind == "p" and y != length:
            raise ValueError("a p-BFB cycle requires an empty R")
        if operation.kind in ("i", "f") and y == length:
            raise ValueError(
                f"a {operation.kind}-BFB cycle requires a nonempty R"
            )
        if operation.kind not in ("p", "i", "f"):
            raise ValueError(f"unknown operation kind {operation.kind!r}")
        return State(
            sequence.bfb_target(x, y, operation.is_forward_deletion),
            self.used,
        )

    @classmethod
    def from_string(cls, text, parse_sequence=Sequence.from_string):
        """Parse a serialized state.

        Raises:
            ValueError: If the text is not a serialized state.
        """
        text = text.strip()
        if not text.startswith("{"):
            raise ValueError(f"invalid state {text!r}: missing source set")
        close = text.find("}")
        if close < 0:
            raise ValueError(f"invalid state {text!r}: unterminated source set")
        used = 0
        indices_text = text[1:close]
        if indices_text:
            for index_text in indices_text.split(","):
                index = int(index_text)
                if index < 0:
                    raise ValueError(f"invalid state {text!r}: negative index")
                if used >> index & 1:
                    raise ValueError(
                        f"invalid state {text!r}: repeated source index"
                    )
                used |= 1 << index
        return cls(parse_sequence(text[close + 1:]), used)


def start_state():
    """Return the state every history starts from."""
    return State(Sequence(), 0)


def end_state(sources, target):
    """Return the state every history reaching ``target`` ends in."""
    return State(target, full_mask(sources))


def validate_sources(sources):
    """Return ``sources`` as a tuple after checking that none is empty.

    Raises:
        ValueError: If a source is empty.
    """
    sources = tuple(sources)
    for index, source in enumerate(sources):
        if len(source) == 0:
            raise ValueError(f"source {index} is empty")
    return sources


def total_source_length(sources):
    """Return the summed length of the sources."""
    return sum(len(source) for source in sources)
