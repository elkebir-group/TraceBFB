#!/usr/bin/env python3
"""Differential test harness for the C++ port of NBFB_fast_dist.

Every instance is solved twice, once by ``bfb_dist.nbfb_fast.NBFB_fast_dist``
and once by the ``nbfb`` binary of this directory. The two must agree on the
total distance and on the number of t-BFB cycles, and the history the C++
reports must replay to the target with the Python sequence operations.

Parts:

* ``random``       simulated histories from ``bfb_dist.simulation``
* ``unreachable``  random targets that mostly cannot be built
* ``castle``       the stored longest walks of the four CASTLE loci
* ``all``          all of the above

Usage:

    PYTHONPATH=../src python3 verify.py --part all
    PYTHONPATH=../src python3 verify.py --part random -n 3000
    PYTHONPATH=../src python3 verify.py --part castle --locus chr11q --limit 2000
"""

import argparse
import json
import math
import os
import random
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

from bfb_dist.history import Operation, start_state  # noqa: E402
from bfb_dist.nbfb_fast import NBFB_fast_dist  # noqa: E402
from bfb_dist.sequence import Sequence  # noqa: E402
from bfb_dist.simulation import (  # noqa: E402
    _set_seed,
    simulate_history,
    simulate_sources,
)

NBFB = os.path.join(HERE, "nbfb")
CASTLE = "/home/melkebir/projects/BFB-dist/CASTLE"
WALKS = os.path.join(CASTLE, "derived", "longest_walks")

LOCI = {
    "chr11q": ("HCC1954_chr11q.json", "HCC1954_chr11q.json", [("T", "E")]),
    "chr6q": ("HCC1395_chr6q.json", "HCC1395_chr6q_open.json", []),
    "chr15q": ("HCC1954_chr15q_Pcen.json", "HCC1954_chr15q_closed.json", []),
    "chr5p": ("H1437_chr5p_Qcen.json", "H1437_chr5p.json", []),
}


# ----------------------------------------------------------------- helpers


def parse_ops(text):
    """Parse the ``operations:`` line of the C++ tool."""
    operations = []
    text = text.strip()
    if not text:
        return operations
    for field in text.split(";"):
        operations.append(Operation.from_fields(field.split(",")))
    return operations


def replay(operations, sources, target):
    """Return whether applying ``operations`` from the start reaches ``target``."""
    state = start_state()
    for operation in operations:
        state = state.apply(operation, sources)
    return state.sequence == target


def python_answer(sources, target):
    """Return ``(t_cycles, total)`` of the Python algorithm, with inf allowed."""
    graph = NBFB_fast_dist(sources, target)
    distance = graph.calculate_bfb_distance()
    if math.isinf(distance):
        return ("inf", "inf")
    return (str(int(graph.translocation_count())), str(int(distance)))


def run_single(path, max_mem_gb):
    """Run the C++ tool on one instance file and return its parsed output."""
    result = subprocess.run(
        [NBFB, "-i", path, "--history", "--max-mem-gb", str(max_mem_gb)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"nbfb failed on {path}: {result.stderr.strip()}")
    lines = result.stdout.rstrip("\n").split("\n")
    return {
        "t_cycles": lines[0].split(": ", 1)[1].strip(),
        "types": lines[1].split(": ", 1)[1].strip(),
        "ops": parse_ops(lines[2].split(": ", 1)[1]),
        "total": lines[-1].strip(),
    }


class Counter:
    def __init__(self):
        self.checked = 0
        self.mismatches = 0
        self.replays = 0
        self.replay_failures = 0
        self.details = []

    def record(self, label, sources, target, expected, got, operations):
        self.checked += 1
        if expected != got:
            self.mismatches += 1
            if len(self.details) < 10:
                self.details.append(
                    f"{label}: sources={[str(s) for s in sources]} "
                    f"target={target} python={expected} cpp={got}"
                )
            return
        if got[1] == "inf":
            return
        self.replays += 1
        if not replay(operations, sources, target):
            self.replay_failures += 1
            if len(self.details) < 10:
                self.details.append(
                    f"{label}: history does not replay to the target: "
                    f"{[str(o) for o in operations]}"
                )

    def report(self, name, extra=""):
        print(
            f"  {name}: checked={self.checked} mismatches={self.mismatches} "
            f"replays_ok={self.replays - self.replay_failures}/{self.replays}"
            f"{extra}"
        )
        for line in self.details:
            print(f"    ! {line}")
        return self.mismatches == 0 and self.replay_failures == 0


def write_instance(path, sources, target):
    with open(path, "w", encoding="utf-8") as handle:
        for source in sources:
            handle.write(f"{source}\n")
        handle.write(f"{target}\n")


# ------------------------------------------------------------ random part


def block_sources(sources):
    """Every contiguous block of every source in both orientations.

    The same source set shape the CASTLE runs use, which lets a small
    instance exercise source reuse and a source set much larger than the
    history that produced the target.
    """
    blocks = []
    seen = set()
    for source in sources:
        for i in range(len(source)):
            for j in range(i + 1, len(source) + 1):
                for block in (source[i:j], source[i:j].reverse_complement()):
                    if block.segments not in seen:
                        seen.add(block.segments)
                        blocks.append(block)
    return blocks


def random_instance(rng, seed):
    """Return one simulated instance, or None when the shape is rejected."""
    n_sources = rng.randint(1, 4)
    n_segments = rng.randint(1, 4)
    total_segments = n_sources * n_segments
    if not 3 <= total_segments <= 12:
        return None
    cycles = rng.randint(n_sources, 8)
    short_deletions = rng.random() < 0.5
    _set_seed(seed)
    sources = simulate_sources(n_sources, n_segments)
    steps, _ = simulate_history(sources, cycles, short_deletions)
    target = steps[-1][0].sequence
    if rng.random() < 0.34:
        # Source reuse: hand the solver every contiguous block instead of
        # the sources the history used.
        return block_sources(sources), target, short_deletions
    return sources, target, short_deletions


def part_random(wanted, max_tokens, max_mem_gb, scratch):
    rng = random.Random(20260919)
    counter = Counter()
    path = os.path.join(scratch, "random.txt")
    generated = 0
    skipped = 0
    short_count = 0
    python_seconds = 0.0
    cpp_seconds = 0.0
    seed = 0
    while counter.checked < wanted:
        seed += 1
        generated += 1
        instance = random_instance(rng, seed)
        if instance is None:
            continue
        sources, target, short_deletions = instance
        if len(target) > max_tokens:
            skipped += 1
            continue
        short_count += int(short_deletions)
        started = time.perf_counter()
        expected = python_answer(sources, target)
        python_seconds += time.perf_counter() - started
        write_instance(path, sources, target)
        started = time.perf_counter()
        output = run_single(path, max_mem_gb)
        cpp_seconds += time.perf_counter() - started
        counter.record(
            f"random seed={seed}",
            sources,
            target,
            expected,
            (output["t_cycles"], output["total"]),
            output["ops"],
        )
    extra = (
        f" (from {generated} draws, {skipped} over {max_tokens} tokens, "
        f"{short_count} with short_deletions; python {python_seconds:.1f}s, "
        f"cpp+spawn {cpp_seconds:.1f}s)"
    )
    return counter.report("random simulated", extra)


def part_reference(wanted, max_mem_gb, scratch):
    """Compare the C++, the fast Python and the O(n^6) reference search.

    The reference ``NBFB_dist`` enumerates strings, so this part only runs on
    tiny instances. A disagreement here would mean the index-tuple rules
    themselves are wrong, in the Python as much as in the port.
    """
    from bfb_dist.nbfb import NBFB_dist  # noqa: E402

    rng = random.Random(4242)
    counter = Counter()
    path = os.path.join(scratch, "reference.txt")
    reference_mismatches = 0
    checked = 0
    seed = 0
    while checked < wanted:
        seed += 1
        n_sources = rng.randint(1, 2)
        cycles = rng.randint(n_sources, n_sources + 3)
        _set_seed(50000 + seed)
        sources = simulate_sources(n_sources, rng.randint(1, 2))
        steps, _ = simulate_history(sources, cycles, rng.random() < 0.5)
        target = steps[-1][0].sequence
        if len(target) > 12:
            continue
        if rng.random() < 0.5:
            sources = block_sources(sources)
        checked += 1
        reference = NBFB_dist(sources, target)
        reference_distance = reference.calculate_bfb_distance()
        if math.isinf(reference_distance):
            expected = ("inf", "inf")
        else:
            expected = (
                str(int(reference.translocation_count())),
                str(int(reference_distance)),
            )
        fast = python_answer(sources, target)
        if fast != expected:
            reference_mismatches += 1
            if len(counter.details) < 10:
                counter.details.append(
                    f"reference seed={seed}: NBFB_dist={expected} "
                    f"NBFB_fast_dist={fast} "
                    f"sources={[str(s) for s in sources]} target={target}"
                )
        write_instance(path, sources, target)
        output = run_single(path, max_mem_gb)
        counter.record(
            f"reference seed={seed}",
            sources,
            target,
            expected,
            (output["t_cycles"], output["total"]),
            output["ops"],
        )
    ok = counter.report(
        "reference NBFB_dist",
        f" (fast-vs-reference mismatches={reference_mismatches})",
    )
    return ok and reference_mismatches == 0


def part_unreachable(wanted, max_mem_gb, scratch):
    rng = random.Random(777)
    counter = Counter()
    path = os.path.join(scratch, "unreachable.txt")
    infinite = 0
    for trial in range(wanted):
        _set_seed(90000 + trial)
        n_sources = rng.randint(1, 3)
        sources = simulate_sources(n_sources, rng.randint(1, 3))
        alphabet = []
        for source in sources:
            for segment in source.segments:
                alphabet.append(segment)
                alphabet.append(segment.reverse_complement())
        if rng.random() < 0.25:
            # A target using a segment that occurs in no source.
            stranger = simulate_sources(1, 1)[0]
            alphabet.append(stranger[0])
        length = rng.randint(2, 12)
        target = Sequence(rng.choice(alphabet) for _ in range(length))
        expected = python_answer(sources, target)
        infinite += int(expected[1] == "inf")
        write_instance(path, sources, target)
        output = run_single(path, max_mem_gb)
        counter.record(
            f"unreachable trial={trial}",
            sources,
            target,
            expected,
            (output["t_cycles"], output["total"]),
            output["ops"],
        )
    return counter.report(
        "random targets", f" ({infinite} of {counter.checked} unreachable)"
    )


# ------------------------------------------------------------ castle part


def castle_instance(model_name, appends):
    """Rebuild the source set and coordinates of longest_walks.py (_init/_tok)."""
    sys.path.insert(0, os.path.join(CASTLE, "scripts"))
    from molecules import build_chains  # noqa: E402

    model = json.load(
        open(os.path.join(CASTLE, "models", model_name), encoding="utf-8")
    )
    chains = build_chains(model, appends)
    coords = {}
    for chain_index, chain in enumerate(chains):
        for position, name in enumerate(chain):
            coords[name.upper()] = (
                f"chr{chain_index}",
                1000 * (position + 1),
                1000 * (position + 1) + 900,
            )

    def token(name):
        chromosome, low, high = coords[name.upper()]
        if name[0].isupper():
            return f"{chromosome}:{low}-{high}[R]"
        return f"{chromosome}:{high}-{low}[C]"

    blocks = []
    for chain in chains:
        for i in range(len(chain)):
            for j in range(i + 1, len(chain) + 1):
                block = list(chain[i:j])
                blocks.append(block)
                blocks.append([x.swapcase() for x in reversed(block)])
    sources = [
        Sequence.from_string(";".join(token(c) for c in block))
        for block in blocks
    ]
    return sources, token


def part_castle(locus, limit, threads, max_mem_gb, scratch, time_python):
    walk_file, model_name, appends = LOCI[locus]
    data = json.load(open(os.path.join(WALKS, walk_file), encoding="utf-8"))
    results = data["results"]
    if limit and len(results) > limit:
        results = results[:limit]
    sources, token = castle_instance(model_name, appends)

    sources_path = os.path.join(scratch, f"{locus}.sources.txt")
    targets_path = os.path.join(scratch, f"{locus}.targets.txt")
    out_path = os.path.join(scratch, f"{locus}.out.tsv")
    with open(sources_path, "w", encoding="utf-8") as handle:
        for source in sources:
            handle.write(f"{source}\n")
    targets = []
    with open(targets_path, "w", encoding="utf-8") as handle:
        for record in results:
            letters = record["walk"].split(" ")
            target = Sequence.from_string(";".join(token(c) for c in letters))
            targets.append(target)
            handle.write(f"{target}\n")

    started = time.perf_counter()
    process = subprocess.run(
        [
            NBFB,
            "--sources",
            sources_path,
            "--targets",
            targets_path,
            "-j",
            str(threads),
            "-o",
            out_path,
            "--history",
            "--max-mem-gb",
            str(max_mem_gb),
        ],
        capture_output=True,
        text=True,
    )
    wall = time.perf_counter() - started
    if process.returncode != 0:
        print(f"  {locus}: nbfb failed: {process.stderr.strip()}")
        return False

    counter = Counter()
    cpu_seconds = 0.0
    with open(out_path, encoding="utf-8") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle][1:]
    if len(rows) != len(results):
        print(f"  {locus}: {len(rows)} rows for {len(results)} targets")
        return False
    for record, row, target in zip(results, rows, targets):
        stored = record.get("nbfb")
        types = record.get("nbfb_types", "")
        if not isinstance(stored, int):
            continue
        expected = (str(types.count("t")), str(stored))
        cpu_seconds += float(row[5])
        counter.record(
            f"{locus} walk={record['walk'][:40]}",
            sources,
            target,
            expected,
            (row[2], row[3]),
            parse_ops(row[6]) if len(row) > 6 else [],
        )

    python_stored = sum(record.get("nbfb_s", 0.0) or 0.0 for record in results)
    extra = (
        f" (wall {wall:.1f}s at -j {threads}, cpp cpu {cpu_seconds:.1f}s, "
        f"stored python {python_stored:.1f}s)"
    )
    ok = counter.report(f"CASTLE {locus}", extra)

    if time_python:
        # Re-time Python on as many walks as fit in about time_python
        # seconds, and time the same walks with one C++ worker, so that the
        # comparison is single core on both sides.
        started = time.perf_counter()
        used = 0
        for target in targets:
            NBFB_fast_dist(sources, target)
            used += 1
            if time.perf_counter() - started >= time_python:
                break
        python_seconds = time.perf_counter() - started
        sample_path = os.path.join(scratch, f"{locus}.sample.txt")
        with open(sample_path, "w", encoding="utf-8") as handle:
            for target in targets[:used]:
                handle.write(f"{target}\n")
        started = time.perf_counter()
        subprocess.run(
            [NBFB, "--sources", sources_path, "--targets", sample_path,
             "-j", "1", "-o", os.devnull, "--max-mem-gb", str(max_mem_gb)],
            check=True,
            capture_output=True,
        )
        cpp_seconds = time.perf_counter() - started
        print(
            f"    single-core timing over {used} walks: "
            f"python {python_seconds:.2f}s "
            f"({1000 * python_seconds / used:.1f} ms/walk), "
            f"cpp {cpp_seconds:.2f}s "
            f"({1000 * cpp_seconds / used:.1f} ms/walk), "
            f"speedup {python_seconds / max(cpp_seconds, 1e-9):.0f}x"
        )
    return ok


# ------------------------------------------------------------------- main


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--part",
        default="all",
        choices=["all", "random", "unreachable", "reference", "castle"],
    )
    parser.add_argument("-n", "--count", type=int, default=3000)
    parser.add_argument("--unreachable-count", type=int, default=300)
    parser.add_argument("--reference-count", type=int, default=200)
    parser.add_argument("--max-tokens", type=int, default=120)
    parser.add_argument("--locus", default="")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("-j", "--threads", type=int, default=1)
    parser.add_argument("--max-mem-gb", type=float, default=2.0)
    parser.add_argument(
        "--time-python",
        type=int,
        default=0,
        help="seconds of Python time to spend per locus on a single-core "
        "timing comparison (0 disables it)",
    )
    parser.add_argument(
        "--scratch",
        default=os.environ.get("NBFB_SCRATCH", "/tmp"),
    )
    args = parser.parse_args()

    if not os.path.exists(NBFB):
        sys.exit(f"build the tool first: make -C {HERE}")
    os.makedirs(args.scratch, exist_ok=True)

    ok = True
    print("verification")
    if args.part in ("all", "random"):
        ok &= part_random(
            args.count, args.max_tokens, args.max_mem_gb, args.scratch
        )
    if args.part in ("all", "unreachable"):
        ok &= part_unreachable(
            args.unreachable_count, args.max_mem_gb, args.scratch
        )
    if args.part in ("all", "reference"):
        ok &= part_reference(
            args.reference_count, args.max_mem_gb, args.scratch
        )
    if args.part in ("all", "castle"):
        loci = [args.locus] if args.locus else list(LOCI)
        defaults = {
            "chr11q": 2000,
            "chr6q": 200,
            "chr15q": 200,
            "chr5p": 30,
        }
        for locus in loci:
            ok &= part_castle(
                locus,
                args.limit or defaults[locus],
                args.threads,
                args.max_mem_gb,
                args.scratch,
                args.time_python,
            )
    print("ALL CHECKS PASSED" if ok else "FAILURES")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
