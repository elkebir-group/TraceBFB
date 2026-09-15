"""Random source and BFB-history simulation, including translocations."""

import random

from bfb_dist.history import Operation, start_state
from bfb_dist.sequence import Segment, Sequence
from bfb_dist.utils import CHROMOSOME_LENGTHS, format_operation_line


def _set_seed(seed):
    """Seed the pseudo-random generator used by simulation functions."""
    random.seed(seed)


def simulate_reference_sequence(
    n_segments,
    min_segment_len=1000,
    max_segment_len=5000,
    starting_position=None,
    chromosome=None,
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

    if chromosome is None:
        eligible_chromosomes = [
            name
            for name, length in CHROMOSOME_LENGTHS.items()
            if length >= minimum_chromosome_length
        ]
        if not eligible_chromosomes:
            raise ValueError(
                "no chromosome can contain the requested reference sequence"
            )
        chromosome = random.choice(eligible_chromosomes)
    elif CHROMOSOME_LENGTHS.get(chromosome, 0) < minimum_chromosome_length:
        raise ValueError(
            f"{chromosome} cannot contain the requested reference sequence"
        )

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


def simulate_sources(
    n_sources,
    n_segments,
    min_segment_len=1000,
    max_segment_len=5000,
):
    """Generate ``n_sources`` reference sequences on distinct chromosomes.

    Distinct chromosomes make every source a translocation partner of the
    others. When more sources than eligible chromosomes are requested,
    chromosomes are reused.

    Raises:
        ValueError: If ``n_sources`` is not positive or the sequence
            arguments are invalid.
    """
    if n_sources <= 0:
        raise ValueError("n_sources must be positive")
    required_span = n_segments * max_segment_len
    eligible = [
        name
        for name, length in CHROMOSOME_LENGTHS.items()
        if length >= required_span
    ]
    if not eligible:
        raise ValueError(
            "no chromosome can contain the requested reference sequence"
        )
    if n_sources <= len(eligible):
        chromosomes = random.sample(eligible, n_sources)
    else:
        chromosomes = [random.choice(eligible) for _ in range(n_sources)]
    return [
        simulate_reference_sequence(
            n_segments, min_segment_len, max_segment_len, chromosome=name
        )
        for name in chromosomes
    ]


def _random_fold_back(sequence, max_overhang=None):
    """Return a random p-, i-, or f-BFB operation on ``sequence``.

    When ``max_overhang`` is given, an f-BFB cycle removes fewer than
    ``max_overhang`` segments, so that the history satisfies the short
    deletions hypothesis of the paper and is therefore nested.
    """
    length = len(sequence)
    x = random.randint(0, length - 1)
    y = random.randint(x + 1, length)
    if y == length:
        return Operation("p", x, y)
    if random.random() <= 0.5:
        return Operation("i", x, y)
    if max_overhang is not None and length - y >= max_overhang:
        if max_overhang <= 1:
            return Operation("i", x, y)
        y = random.randint(max(x + 1, length - max_overhang + 1), length)
        if y == length:
            return Operation("p", x, y)
    return Operation("f", x, y)


def simulate_history(sources, num_cycles, short_deletions=False, color_stream=False):
    """Simulate a random history over ``sources`` of ``num_cycles`` cycles.

    The history starts from the empty sequence. Its first cycle is a t-BFB
    cycle, and every source is appended by exactly one t-BFB cycle at a
    random position; the remaining cycles are random p-, i-, or f-BFB
    cycles. When ``short_deletions`` is true, every f-BFB cycle removes
    fewer segments than the preceding cycle added, so the history is
    nested.

    Returns:
        A pair ``(steps, ansi_output)``. ``steps`` is the list of
        ``(state, operation)`` pairs followed by the final state paired
        with ``None``. ``ansi_output`` is a colored sequence-per-line
        representation (without the empty start) when ``color_stream`` is
        true, and otherwise ``None``.

    Raises:
        ValueError: If there is no source, or fewer cycles than sources
            are requested.
    """
    sources = list(sources)
    if not sources:
        raise ValueError("at least one source is required")
    if any(len(source) == 0 for source in sources):
        raise ValueError("sources cannot be empty")
    if num_cycles < len(sources):
        raise ValueError("num_cycles cannot be smaller than the number of sources")

    order = list(range(len(sources)))
    random.shuffle(order)
    positions = {0: order[0]}
    if len(sources) > 1:
        later = random.sample(range(1, num_cycles), len(sources) - 1)
        for position, index in zip(sorted(later), order[1:]):
            positions[position] = index

    state = start_state()
    steps = []
    previous_gain = None
    for cycle in range(num_cycles):
        if cycle in positions:
            operation = Operation.translocation(positions[cycle])
            gain = len(sources[positions[cycle]])
        else:
            operation = _random_fold_back(
                state.sequence,
                previous_gain if short_deletions else None,
            )
            gain = operation.y - operation.x
        steps.append((state, operation))
        state = state.apply(operation, sources)
        previous_gain = gain
    steps.append((state, None))

    color_output = None
    if color_stream:
        lines = [
            format_operation_line(step_state.sequence, operation, colored=True)
            for step_state, operation in steps[1:-1]
        ]
        lines.append(str(state.sequence))
        color_output = "\n".join(lines) + "\n"

    return steps, color_output
