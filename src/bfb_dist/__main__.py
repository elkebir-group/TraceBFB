import argparse
import sys
from bfb_dist.sequence import Sequence
from bfb_dist.simulation import simulate_reference_sequence, simulate_bfb
from bfb_dist.distance import BFBDist

def handle_simulate(args):
    ref = simulate_reference_sequence(args.n_segments, args.min_len, args.max_len)
    _, _, ops = simulate_bfb(ref, args.cycles, args.verbose)
    with open(args.output, 'w') as f:
        for seq in ops:
            f.write(f"{seq}\n")
    print(f"Simulation saved to {args.output}")

def handle_bfbdist(args):
    try:
        with open(args.input, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        sys.exit(f"Error: Input file '{args.input}' not found.")

    if len(lines) < 2:
        sys.exit("Error: Input file must contain at least 2 lines (source and target).")

    source = Sequence.from_string(lines[0])
    target = Sequence.from_string(lines[1])
    G, _, min_cost, paths = BFBDist(source, target)
    print(f"Minimum Cost: {min_cost}")

    with open(args.output, 'w') as f:
        f.write("=== GRAPH EDGES ===\n")
        for parent, children in G.items():
            for child in children:
                f.write(f"{parent} -> {child}\n")
        f.write("\n=== PATHS ===\n")
        for idx, path in enumerate(paths, 1):
            f.write(f"Path {idx}: {' -> '.join(map(str, path))}\n")
    print(f"Graph and paths saved to {args.output}")

def main():
    parser = argparse.ArgumentParser(description="BFB-dist command line utility")
    subparsers = parser.add_subparsers(dest="command", required=True)

    sim_parser = subparsers.add_parser("simulate")
    sim_parser.add_argument("-n", "--n-segments", type=int, default=5)
    sim_parser.add_argument("--min-len", type=int, default=1000)
    sim_parser.add_argument("--max-len", type=int, default=5000)
    sim_parser.add_argument("-c", "--cycles", type=int, default=3)
    sim_parser.add_argument("-o", "--output", type=str, required=True)
    sim_parser.add_argument("-v", "--verbose", action="store_true")
    sim_parser.set_defaults(func=handle_simulate)

    dist_parser = subparsers.add_parser("bfbdist")
    dist_parser.add_argument("-i", "--input", type=str, required=True)
    dist_parser.add_argument("-o", "--output", type=str, required=True)
    dist_parser.set_defaults(func=handle_bfbdist)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()