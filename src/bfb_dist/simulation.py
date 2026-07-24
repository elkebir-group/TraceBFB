"""Random reference-sequence and BFB-cycle simulation."""

import random

from bfb_dist.sequence import Segment, Sequence
from bfb_dist.utils import CHROMOSOME_LENGTHS, format_sequence


def _set_seed(seed):
    """Seed the pseudo-random generator used by simulation functions."""
    random.seed(seed)


def simulate_reference_sequence(
    n_segments,
    min_segment_len=1000,
    max_segment_len=5000,
    starting_position=None,
):
    """Generate a contiguous reference-oriented genomic sequence.

    A chromosome and, unless provided, a starting position are selected
    randomly. Segment lengths are sampled inclusively between the requested
    bounds.

    Raises:
        ValueError: If the arguments are invalid or no chromosome can contain
            the requested maximum span.
    """
    if n_segments <= 0:
        raise ValueError("n_segments must be positive")
    if min_segment_len <= 0 or max_segment_len <= 0:
        raise ValueError("segment lengths must be positive")
    if min_segment_len > max_segment_len:
        raise ValueError(
            "min_segment_len cannot be greater than max_segment_len"
        )
    if starting_position is not None and starting_position < 0:
        raise ValueError("starting_position cannot be negative")

    required_span = n_segments * max_segment_len
    minimum_chromosome_length = required_span
    if starting_position is not None:
        minimum_chromosome_length += starting_position

    eligible_chromosomes = [
        chromosome
        for chromosome, length in CHROMOSOME_LENGTHS.items()
        if length >= minimum_chromosome_length
    ]
    if not eligible_chromosomes:
        raise ValueError(
            "no chromosome can contain the requested reference sequence"
        )

    chromosome = random.choice(eligible_chromosomes)
    chromosome_length = CHROMOSOME_LENGTHS[chromosome]
    x = (
        random.randint(0, chromosome_length - required_span)
        if starting_position is None
        else starting_position
    )

    ref = Sequence()
    for _ in range(n_segments):
        y = x + random.randint(min_segment_len, max_segment_len)
        ref += Segment(chromosome, x, y, False)
        x = y
    return ref


def simulate_bfb(reference, num_cycles, color_stream=False):
    """Simulate successive random BFB operations.

    Returns:
        A pair ``(sequences, ansi_output)``. ``sequences`` contains the
        reference followed by every generated target. ``ansi_output`` is a
        sequence-per-line colored representation when ``color_stream`` is
        true, and otherwise is ``None``.

    Raises:
        ValueError: If the cycle count is negative or a nonzero number of
            cycles is requested for an empty reference.
    """
    if num_cycles < 0:
        raise ValueError("num_cycles cannot be negative")
    if num_cycles > 0 and not reference:
        raise ValueError("reference cannot be empty when simulating BFB cycles")

    template = reference
    ops = [reference]
    color_lines = [] if color_stream else None

    for cycle in range(num_cycles):
        source = template
        x = random.randint(0, len(source) - 1)
        y = random.randint(x + 1, len(source))
        is_forward_deletion = random.random() <= 0.5

        if color_stream:
            decomposition = source.decompose(x, y)
            color_lines.append(
                format_sequence(decomposition, colored=True)
            )

        template = source.bfb_target(x, y, is_forward_deletion)
        ops.append(template)

    color_output = None
    if color_stream:
        color_lines.append(str(template))
        color_output = "\n".join(color_lines) + "\n"

    return ops, color_output
