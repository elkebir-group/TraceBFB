"""Command-line interface for BFBdist."""

import argparse
import sys


def handle_simulate(args):
    """Run simulation and write its requested path artifacts."""
    from bfb_dist.simulation import (
        _set_seed,
        simulate_history,
        simulate_sources,
    )

    if args.seed != -1:
        _set_seed(args.seed)

    try:
        sources = simulate_sources(
            args.n_sources,
            args.n_segments,
            args.min_len,
            args.max_len,
        )
        steps, color_output = simulate_history(
            sources,
            args.cycles,
            args.short_deletions,
            args.colored,
        )
    except ValueError as error:
        sys.exit(f"Error: {error}")

    target = steps[-1][0].sequence
    try:
        with open(f"{args.out}.input", "w", encoding="utf-8") as output:
            for source in sources:
                output.write(f"{source}\n")
            output.write(f"{target}\n")
        with open(f"{args.out}.paths", "w", encoding="utf-8") as output:
            for state, _ in steps[1:]:
                output.write(f"{state.sequence}\n")
        if args.colored:
            with open(
                f"{args.out}.paths.ansi",
                "w",
                encoding="utf-8",
            ) as output:
                output.write(color_output)
    except OSError as error:
        sys.exit(f"Error writing output file: {error}")


def read_instance(path):
    """Read sources and target from ``path``.

    Every line but the last is a source and the last line is the target.

    Raises:
        ValueError: If the file has fewer than two lines, contains a blank
            line, or contains a malformed sequence.
    """
    from bfb_dist.sequence import Sequence

    with open(path, "r", encoding="utf-8") as input_file:
        lines = input_file.read().splitlines()

    if len(lines) < 2 or any(not line.strip() for line in lines):
        raise ValueError(
            "input file must contain one or more sources followed by the "
            "target, one nonempty sequence per line"
        )
    sequences = [Sequence.from_string(line) for line in lines]
    return sequences[:-1], sequences[-1]


def handle_bfbdist(args):
    """Infer a BFB graph and write graph, DOT, and path artifacts."""
    from bfb_dist.distance import BFB_dist
    from bfb_dist.nbfb import NBFB_dist

    try:
        sources, target = read_instance(args.input)
    except OSError as error:
        sys.exit(f"Error reading input file: {error}")
    except ValueError as error:
        sys.exit(f"Error: {error}")

    compute = {
        "bfb": BFB_dist,
        "nbfb": NBFB_dist,
    }[args.mode]
    graph = compute(sources, target)

    try:
        with open(f"{args.out}.graph", "w", encoding="utf-8") as output:
            graph.write_text(output)
        with open(f"{args.out}.dot", "w", encoding="utf-8") as output:
            graph.write_dot(output)
        with open(f"{args.out}.paths", "w", encoding="utf-8") as output:
            output.write(graph.write_paths(colored=False))
        if args.colored:
            with open(
                f"{args.out}.paths.ansi",
                "w",
                encoding="utf-8",
            ) as output:
                output.write(graph.write_paths(colored=True))
    except OSError as error:
        sys.exit(f"Error writing output file: {error}")

    print(graph.calculate_bfb_distance())


def main():
    """Parse command-line arguments and dispatch the selected command."""
    parser = argparse.ArgumentParser(description="BFBdist command line utility")
    subparsers = parser.add_subparsers(dest="command", required=True)

    sim_parser = subparsers.add_parser("simulate")
    sim_parser.add_argument("-s", "--seed", type=int, default=-1)
    sim_parser.add_argument("-m", "--n-sources", type=int, default=1)
    sim_parser.add_argument("-n", "--n-segments", type=int, default=5)
    sim_parser.add_argument("--min-len", type=int, default=1000)
    sim_parser.add_argument("--max-len", type=int, default=5000)
    sim_parser.add_argument("-c", "--cycles", type=int, default=3)
    sim_parser.add_argument("--short-deletions", action="store_true")
    sim_parser.add_argument("-o", "--out", type=str, required=True)
    sim_parser.add_argument("--colored", action="store_true")
    sim_parser.set_defaults(func=handle_simulate)

    dist_parser = subparsers.add_parser("distance")
    dist_parser.add_argument("-i", "--input", type=str, required=True)
    dist_parser.add_argument("-o", "--out", type=str, required=True)
    dist_parser.add_argument("--colored", action="store_true")
    dist_parser.add_argument(
        "--mode",
        choices=["bfb", "nbfb"],
        default="bfb",
    )
    dist_parser.set_defaults(func=handle_bfbdist)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
