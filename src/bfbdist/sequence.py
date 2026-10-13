"""Sequence and segment representations used by BFBdist."""

import re
from dataclasses import dataclass


_SEGMENT_PATTERN = re.compile(
    r"(?P<chromosome>[^:;\s]+):(?P<start>\d+)-(?P<end>\d+)"
    r"\[(?P<strand>[RC])\]"
)

_SENTINEL = object()


def _z_function(elements):
    """Return the Z-array of ``elements`` (Gusfield's Z-algorithm).

    ``z[i]`` is the length of the longest run starting at ``elements[i]``
    that matches the run starting at ``elements[0]``. Linear time, over any
    sequence of elements supporting ``==``, not just characters.
    """
    length = len(elements)
    z = [0] * length
    left = right = 0
    for i in range(1, length):
        if i < right:
            z[i] = min(right - i, z[i - left])
        while (
            i + z[i] < length
            and elements[z[i]] == elements[i + z[i]]
        ):
            z[i] += 1
        if i + z[i] > right:
            left, right = i, i + z[i]
    return z


@dataclass(frozen=True)
class Segment:
    """An immutable oriented genomic segment."""

    chr: str
    start: int
    end: int
    complement: bool

    def __str__(self):
        """Serialize the segment as chromosome, coordinates, and strand."""
        strand = "C" if self.complement else "R"
        return f"{self.chr}:{self.start}-{self.end}[{strand}]"

    __repr__ = __str__

    def reverse_complement(self):
        """Return the segment in the opposite orientation."""
        return Segment(self.chr, self.end, self.start, not self.complement)


@dataclass(frozen=True)
class LetterAsSegment:
    """An immutable letter used as a segment in symbolic test sequences."""

    letter: str

    def __str__(self):
        """Return the represented letter."""
        return self.letter

    __repr__ = __str__

    def reverse_complement(self):
        """Return the opposite-case representation of the letter."""
        return LetterAsSegment(self.letter.swapcase())


@dataclass(frozen=True)
class Sequence:
    """An immutable ordered collection of oriented segments."""

    segments: tuple = ()

    def __init__(self, segments=None):
        """Create a sequence from an iterable of segments."""
        object.__setattr__(
            self,
            "segments",
            tuple(segments) if segments is not None else (),
        )

    def __str__(self):
        """Serialize the sequence as semicolon-separated segments."""
        return ";".join(map(str, self.segments))

    __repr__ = __str__

    def __getitem__(self, index):
        """Return one segment, or a Sequence when given a slice."""
        if isinstance(index, slice):
            return Sequence(self.segments[index])
        return self.segments[index]

    def __len__(self):
        """Return the number of segments."""
        return len(self.segments)

    def __add__(self, other):
        """Concatenate another sequence or append one segment."""
        if isinstance(other, Sequence):
            return Sequence(self.segments + other.segments)
        if isinstance(other, (Segment, LetterAsSegment)):
            return Sequence(self.segments + (other,))
        return NotImplemented

    def reverse_complement(self):
        """Return the reverse complement of the entire sequence."""
        return Sequence(
            segment.reverse_complement()
            for segment in reversed(self.segments)
        )

    def is_prefix_of(self, superseq):
        """Return whether this sequence is a prefix of ``superseq``."""
        return self.segments == superseq.segments[:len(self.segments)]

    def is_suffix_of(self, superseq):
        """Return whether this sequence is a suffix of ``superseq``."""
        length = len(self.segments)
        if length == 0:
            return True
        return self.segments == superseq.segments[-length:]

    def strip_suffix(self, suffix):
        """Return the sequence without ``suffix``, or ``None`` if absent."""
        if not suffix.is_suffix_of(self):
            return None
        return self[:len(self) - len(suffix)]

    def is_infix_of(self, other):
        """Return whether this sequence occurs as a contiguous run of ``other``.

        Linear-time Z-algorithm search over segment tuples (not just
        characters), matching the complexity Theorem thm:star-runtime
        assumes for the foldability test. ``_SENTINEL`` is a fresh object
        that no ``Segment``/``LetterAsSegment`` can ever equal, so it safely
        separates the two runs without needing an alphabet-specific marker.
        """
        n, m = len(self), len(other)
        if n == 0:
            return True
        if n > m:
            return False

        combined = self.segments + (_SENTINEL,) + other.segments
        z = _z_function(combined)
        return any(value >= n for value in z[n + 1:])

    def lcp_length(self, other):
        """Return the length of the longest common prefix with ``other``."""
        length = 0
        for a, b in zip(self.segments, other.segments):
            if a != b:
                break
            length += 1
        return length

    def decompose(self, x, y):
        """Return P, Q, and R where ``self = P + Q + R``."""
        return self[:x], self[x:y], self[y:]

    def bfb_target(self, x, y, is_forward_deletion):
        """Return the target of the specified BFB operation."""
        p, q, r = self.decompose(x, y)
        if is_forward_deletion:
            return (
                p
                + q
                + r.reverse_complement()
                + q.reverse_complement()
            )
        return p + q + r + q.reverse_complement()

    def previous_bfb(self):
        """Return all predecessor sequences and operations producing this one.

        Each returned tuple contains ``(source, x, y,
        is_forward_deletion)``.
        """
        predecessors = set()
        target_length = len(self)

        for source_length in range(
            (target_length + 1) // 2,
            target_length,
        ):
            q = self[source_length:].reverse_complement()
            q_length = target_length - source_length

            for x in range(0, 2 * source_length - target_length + 1):
                y = x + q_length
                if self[x:y] != q:
                    continue

                source = self[:source_length]
                predecessors.add((source, x, y, False))

                forward_deletion_source = (
                    self[:y]
                    + self[y:source_length].reverse_complement()
                )
                predecessors.add(
                    (forward_deletion_source, x, y, True)
                )

        return predecessors

    def successors(self):
        """Return all successor sequences and operations producible from this one.

        Each returned tuple contains ``(child, x, y, is_forward_deletion)``.
        """
        successors = set()
        length = len(self)
        for x in range(length):
            for y in range(x + 1, length + 1):
                successors.add((self.bfb_target(x, y, False), x, y, False))
                if y < length:
                    successors.add((self.bfb_target(x, y, True), x, y, True))
        return successors

    @classmethod
    def from_string(cls, seq_str):
        """Parse a semicolon-separated genomic sequence.

        Raises:
            ValueError: If any serialized segment is malformed.
        """
        seq_str = seq_str.strip()
        if not seq_str:
            return cls()

        segments = []
        for index, serialized_segment in enumerate(
            seq_str.split(";"),
            start=1,
        ):
            match = _SEGMENT_PATTERN.fullmatch(serialized_segment.strip())
            if match is None:
                raise ValueError(
                    "Invalid segment "
                    f"{index}: {serialized_segment!r}"
                )

            segments.append(
                Segment(
                    match.group("chromosome"),
                    int(match.group("start")),
                    int(match.group("end")),
                    match.group("strand") == "C",
                )
            )

        return cls(segments)
