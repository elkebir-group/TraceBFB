"""Command-line interface for TraceBFB."""

import argparse
import sys


def handle_simulate(args):
    """Run simulation and write its requested path artifacts."""
    from tracebfb.simulation import (
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
    from tracebfb.sequence import Sequence

    with open(path, "r", encoding="utf-8") as input_file:
        lines = input_file.read().splitlines()

    if len(lines) < 2 or any(not line.strip() for line in lines):
        raise ValueError(
            "input file must contain one or more sources followed by the "
            "target, one nonempty sequence per line"
        )
    sequences = [Sequence.from_string(line) for line in lines]
    return sequences[:-1], sequences[-1]


def handle_tracebfb(args):
    """Infer a BFB graph and write graph, DOT, and path artifacts."""
    from tracebfb.distance import BFB_dist

    try:
        sources, target = read_instance(args.input)
    except OSError as error:
        sys.exit(f"Error reading input file: {error}")
    except ValueError as error:
        sys.exit(f"Error: {error}")

    graph = BFB_dist(sources, target)

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


def handle_graph(args):
    """Build the breakpoint graph of the regions and write graph.tsv (and subpaths.tsv)."""
    import os

    from tracebfb.graph_construction import build_graph

    reads = []
    for spec in args.reads:
        name, sep, bam = spec.partition("=")
        if not sep or not name or not bam:
            sys.exit(f"Error: --reads expects SET=BAM, got {spec!r}")
        reads.append((name, bam))
    header = {"locus": args.name} if args.name else {}
    try:
        graph, rows = build_graph(
            args.region, args.junctions, args.normal_junctions, args.bam, args.normal_bam,
            args.haploid_depth, cn_files=args.cn, reads=reads, use_subpaths=args.use_subpaths,
            cytoband=args.cytoband, header=header)
    except (OSError, ValueError) as error:
        sys.exit(f"Error: {error}")
    os.makedirs(args.out, exist_ok=True)
    graph.write(os.path.join(args.out, "graph.tsv"))
    if reads:
        with open(os.path.join(args.out, "subpaths.tsv"), "w", encoding="utf-8") as output:
            output.write("string\treads\tread_set\n")
            output.write("".join(f"{w}\t{n}\t{r}\n" for w, n, r in rows))


def main():
    """Parse command-line arguments and dispatch the selected command."""
    parser = argparse.ArgumentParser(description="TraceBFB command line utility")
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
    dist_parser.set_defaults(func=handle_tracebfb)

    graph_parser = subparsers.add_parser(
        "graph", help="breakpoint graph of genomic regions from long reads")
    graph_parser.add_argument("--region", action="append", required=True,
                              help="chrom:start-end, once per window")
    graph_parser.add_argument("--junctions", required=True, help="tumor junction TSV")
    graph_parser.add_argument("--normal-junctions", required=True, help="normal junction TSV")
    graph_parser.add_argument("--bam", required=True, help="tumor BAM (depth)")
    graph_parser.add_argument("--normal-bam", required=True, help="normal BAM (depth)")
    graph_parser.add_argument("--haploid-depth", type=float, required=True,
                              help="tumor read depth of one copy")
    graph_parser.add_argument("--cn", action="append", default=[],
                              help="copy-number BED (copy number in column 5); summed")
    graph_parser.add_argument("--reads", action="append", default=[], metavar="SET=BAM",
                              help="read set for subpath constraints")
    graph_parser.add_argument("--use-subpaths", action="store_true",
                              help="use the subpath constraints as bounds in the fit")
    graph_parser.add_argument("--cytoband", help="UCSC cytoband file (centromere filter)")
    graph_parser.add_argument("--name", help="locus name for the header")
    graph_parser.add_argument("-o", "--out", required=True, help="output directory")
    graph_parser.set_defaults(func=handle_graph)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
