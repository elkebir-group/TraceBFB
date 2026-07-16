class Segment:
    def __init__(self, chromosome, start, end, complement):
        self.chr = chromosome
        self.start = start
        self.end = end
        self.complement = complement

    def __eq__(self, other):
        return isinstance(other, Segment) and (self.chr, self.start, self.end, self.complement) == (other.chr, other.start, other.end, other.complement)

    def __str__(self):
        return f"{self.chr}:{self.start}-{self.end}[{'C' if self.complement else 'R'}]"

    __repr__ = __str__

    def reverse_complement(self):
        return Segment(self.chr, self.end, self.start, not self.complement)


class LetterAsSegment:
    def __init__(self, letter):
        self.letter = letter

    def reverse_complement(self):
        return self.letter.swapcase()


class Sequence:
    def __init__(self, segments=None):
        self.segments = segments if segments is not None else []

    def __eq__(self, other):
        return isinstance(other, Sequence) and self.segments == other.segments

    def __str__(self):
        return ';'.join(map(str, self.segments))

    def __getitem__(self, i):
        return Sequence(self.segments[i]) if isinstance(i, slice) else self.segments[i]

    def __len__(self):
        return len(self.segments)

    def __add__(self, b):
        if isinstance(b, Sequence):
            return Sequence(self.segments + b.segments)
        if isinstance(b, Segment):
            return Sequence(self.segments + [b])
        return NotImplemented

    __repr__ = __str__

    def __hash__(self):
        return hash(str(self))

    def reverse_complement(self):
        return Sequence([s.reverse_complement() for s in reversed(self.segments)])

    def is_prefix_of(self, superseq):
        return self.segments == superseq.segments[:len(self.segments)]

    def previous_bfb(self):
        pot_sources = {}
        l = len(self.segments)
        for i in range((l + 1) // 2, l):
            rc_target = self[i:].reverse_complement()
            target_len = l - i
            for x in range(0, 2 * i - l + 1):
                y = x + target_len
                if self[x:y] == rc_target:
                    pot_sources[self[:i]] = (self[:x], self[x:y], self[y:i])
                    pot_sources[self[:y] + self[y:i].reverse_complement()] = (self[:x], self[x:y], self[y:i].reverse_complement())
        return pot_sources

    @classmethod
    def from_string(cls, seq_str):
        seq_str = seq_str.strip()
        if not seq_str:
            return cls([])
        segments = []
        for seg in seq_str.split(';'):
            chrom, rest = seg.split(':')
            coords, strand_info = rest.split('[')
            start_str, end_str = coords.split('-')
            segments.append(Segment(chrom, int(start_str), int(end_str), strand_info.startswith('C')))
        return cls(segments)