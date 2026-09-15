"""Command-line interface for BFB-dist."""

import argparse
import sys


def handle_simulate(args):
    """Run simulation and write its requested path artifacts."""
    from bfb_dist.simulation import (
        _set_seed,
        simulate_bfb,
        simulate_reference_sequence,
    )

    if args.seed != -1:
        _set_seed(args.seed)

    try:
        ref = simulate_reference_sequence(
            args.n_segments,
            args.min_len,
            args.max_len,
        )
        ops, color_output = simulate_bfb(
            ref,
            args.cycles,
            args.colored,
        )
    except ValueError as error:
        sys.exit(f"Error: {error}")

    try:
        with open(f"{args.out}.paths", "w", encoding="utf-8") as output:
            for sequence in ops:
                output.write(f"{sequence}\n")
        if args.colored:
            with open(
                f"{args.out}.paths.ansi",
                "w",
                encoding="utf-8",
            ) as output:
                output.write(color_output)
    except OSError as error:
        sys.exit(f"Error writing output file: {error}")


def handle_bfbdist(args):
    """Infer a BFB graph and write graph, DOT, and path artifacts."""
    from bfb_dist.distance import BFB_dist
    from bfb_dist.nbfb import NBFB_dist, NBFB_star_dist
    from bfb_dist.sequence import Sequence

    try:
        with open(args.input, "r", encoding="utf-8") as input_file:
            lines = input_file.read().splitlines()
    except OSError as error:
        sys.exit(f"Error reading input file: {error}")

    if len(lines) != 2 or any(not line.strip() for line in lines):
        sys.exit(
            "Error: input file must contain source and target "
            "on two consecutive nonempty lines"
        )

    try:
        source = Sequence.from_string(lines[0])
        target = Sequence.from_string(lines[1])
    except ValueError as error:
        sys.exit(f"Error: {error}")

    compute = {
        "bfb": BFB_dist,
        "nbfb": NBFB_dist,
        "nbfb-star": NBFB_star_dist,
    }[args.mode]
    graph = compute(source, target)

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


def main():
    """Parse command-line arguments and dispatch the selected command."""
    parser = argparse.ArgumentParser(description="BFB-dist command line utility")
    subparsers = parser.add_subparsers(dest="command", required=True)

    sim_parser = subparsers.add_parser("simulate")
    sim_parser.add_argument("-s", "--seed", type=int, default=-1)
    sim_parser.add_argument("-n", "--n-segments", type=int, default=5)
    sim_parser.add_argument("--min-len", type=int, default=1000)
    sim_parser.add_argument("--max-len", type=int, default=5000)
    sim_parser.add_argument("-c", "--cycles", type=int, default=3)
    sim_parser.add_argument("-o", "--out", type=str, required=True)
    sim_parser.add_argument("--colored", action="store_true")
    sim_parser.set_defaults(func=handle_simulate)

    dist_parser = subparsers.add_parser("distance")
    dist_parser.add_argument("-i", "--input", type=str, required=True)
    dist_parser.add_argument("-o", "--out", type=str, required=True)
    dist_parser.add_argument("--colored", action="store_true")
    dist_parser.add_argument(
        "--mode",
        choices=["bfb", "nbfb", "nbfb-star"],
        default="bfb",
    )
    dist_parser.set_defaults(func=handle_bfbdist)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
