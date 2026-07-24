from bfb_dist.sequence import LetterAsSegment, Segment, Sequence


def letters(text):
    return Sequence(LetterAsSegment(letter) for letter in text)


def genomic_sequence(length=3):
    return Sequence(
        Segment("chr1", index * 10, (index + 1) * 10, False)
        for index in range(length)
    )
