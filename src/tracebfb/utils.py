"""Constants and ANSI sequence-formatting helpers."""

CHROMOSOME_LENGTHS = {
    "chr1": 248387328, "chr2": 242696752, "chr3": 201105948, "chr4": 193574945,
    "chr5": 182045439, "chr6": 172126628, "chr7": 160567428, "chr8": 146259331,
    "chr9": 150617247, "chr10": 134758134, "chr11": 135127769, "chr12": 133324548,
    "chr13": 113566686, "chr14": 101161492, "chr15": 99753195, "chr16": 96330374,
    "chr17": 84276897, "chr18": 80542538, "chr19": 61707364, "chr20": 66210255,
    "chr21": 45090682, "chr22": 51324926, "chrX": 154259566, "chrY": 62460029,
    "chrM": 16569,
}


def format_red(s):
    """Wrap a value in bright-red ANSI codes."""
    return f"\033[91m{s}\033[00m"


def format_green(s):
    """Wrap a value in bright-green ANSI codes."""
    return f"\033[92m{s}\033[00m"


def format_yellow(s):
    """Wrap a value in bright-yellow ANSI codes."""
    return f"\033[93m{s}\033[00m"


def format_light_purple(s):
    """Wrap a value in bright-blue ANSI codes."""
    return f"\033[94m{s}\033[00m"


def format_purple(s):
    """Wrap a value in bright-purple ANSI codes."""
    return f"\033[95m{s}\033[00m"


def format_cyan(s):
    """Wrap a value in bright-cyan ANSI codes."""
    return f"\033[96m{s}\033[00m"


def format_light_gray(s):
    """Wrap a value in bright-white ANSI codes."""
    return f"\033[97m{s}\033[00m"


def format_black(s):
    """Wrap a value in bright-black ANSI codes."""
    return f"\033[90m{s}\033[00m"


def format_sequence(parts, colored=False):
    """Format P, Q, and R as one complete sequence.

    When ``colored`` is true, P is red, Q is green, and R is purple.
    """
    p, q, r = parts
    if not colored:
        return str(p + q + r)

    return ";".join(
        formatter(segment)
        for part, formatter in (
            (p, format_red),
            (q, format_green),
            (r, format_purple),
        )
        for segment in part.segments
    )


def format_operation_line(sequence, operation, colored=False):
    """Format the sequence a cycle acts on, colored by the cycle.

    For p-, i-, and f-BFB cycles, P is red, Q is green, and R is purple.
    For a t-BFB cycle, the whole sequence is cyan, as the cycle appends a
    source without decomposing the sequence.
    """
    if not colored:
        return str(sequence)
    if operation.kind == "t":
        return ";".join(format_cyan(segment) for segment in sequence.segments)
    return format_sequence(
        sequence.decompose(operation.x, operation.y),
        colored=True,
    )
