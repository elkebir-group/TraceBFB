"""The breakpoint graph of genomic regions from long reads (`tracebfb graph`).

Inputs: the regions (windows), the tumor and normal junction calls (TSV with the
columns chrA posA sideA chrB posB sideB support; side R when the segment lies to
the left of pos), the tumor and normal BAMs, the haploid depth (tumor read depth
of one copy), optionally copy-number calls (BED, copy number in column 5, summed
over files, e.g. the two haplotypes of Wakhan), and read sets for subpath
constraints.

- Junctions: tumor calls with at least min_support reads, with both ends in the
   regions, not within centromere_margin of a centromere (with a cytoband file),
   not a short deletion (R end then L end, span <= max_small_deletion), and
   somatic (no normal call with both ends within germline_tol, same sides). The
   support of a junction is that of the best-supported tumor call within
   cluster_window of both ends. A junction is a fold-back when its two ends are of
   the same kind on one chromosome within foldback_span, unless an opposite-side
   junction of that shape has both ends within reciprocal_tol (then the two are an
   inversion).
- Segments: every region is cut at the junction ends inside it, and the
   segments get the ids 1, 2, ... in genomic order.
- Depth: bins of the regions; a bin is usable when its normal depth is within
   mapp_tol of the median. A segment with at least 3 usable bins gets the median
   of tumor/normal x median normal; a shorter one the tumor per-base depth,
   normal-corrected when the normal is within 25% of the median.
- Subpath constraints: the oriented segments every read of every read set
   crosses (a read with any alignment in the regions, rebuilt from its SA tag
   from the alignments of MAPQ at least read_mapq; an alignment of at most
   max_gap that overlaps no segment counts as unaligned read), split where a
   block covers no segment, at an unaligned gap of more than max_gap, and at a
   step that is neither a reference edge nor a junction. A constraint is
   written in the reading with more forward segments (ties: the smaller
   reading as integers).
- Fit: integer copy numbers and multiplicities, balanced, every junction at
   least 1, an even number (at least 2) of free ends at chain ends. Objective:
   depth term + Poisson deviance of the junction reads + a prior on the depth
   per copy around the scale of the copy-number calls; the depth per copy is
   scanned over [lo, hi] x the haploid depth. With use_subpaths, every subpath
   constraint is a lower bound. The graph is the fitted one: intact (unrearranged)
   copies of the regions are part of it.
"""

import gzip
import itertools
import math
import random
import re
import statistics
import tempfile
from collections import Counter, defaultdict

import pysam

from tracebfb.graph_file import FOLDBACK_SPAN, GraphFile, GraphSegment, is_foldback, segment_ids

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

DEFAULTS = dict(
    min_support=5, germline_tol=200, max_small_deletion=600, centromere_margin=1_000_000,
    cluster_window=50, foldback_span=FOLDBACK_SPAN, reciprocal_tol=10_000,
    bin=2000, depth_mapq=10, mapp_tol=0.15, trim=2000,
    read_mapq=10, merge=1000, max_gap=2000, pad=1000,
    lo=0.75, hi=1.25, step=0.0125, starts=16, prior_sd=0.10, tol=10_000, seed=1,
)


def chrom_key(c):
    """Genomic order of chromosome names, with or without the chr prefix: 1..22, X, Y, M, rest."""
    n = c[3:] if c.startswith("chr") else c
    if n.isdigit():
        return (0, int(n), "")
    return (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(n, 3), n)


def parse_region(text):
    chrom, span = text.rsplit(":", 1)
    lo, hi = span.replace(",", "").split("-")
    return chrom, int(lo), int(hi)


# ---------------------------------------------------------------- junctions

def read_calls(path):
    """[(end a, end b, support)] with ends (chrom, pos, side), by column name."""
    out, columns = [], None
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            f = raw.rstrip("\n").split("\t")
            if columns is None:
                columns = {name: k for k, name in enumerate(f)}
                missing = {"chrA", "posA", "sideA", "chrB", "posB", "sideB", "support"} - set(columns)
                if missing:
                    raise ValueError(f"{path}: missing columns {sorted(missing)}")
                continue
            get = lambda name: f[columns[name]]  # noqa: E731
            a = (get("chrA"), int(get("posA")), get("sideA"))
            b = (get("chrB"), int(get("posB")), get("sideB"))
            if a[0] == b[0] and a[1] > b[1]:  # the filters read a same-chromosome call left to right
                a, b = b, a
            out.append((a, b, int(get("support"))))
    return out


def centromeres(path):
    cen = {}
    with gzip.open(path, "rt") as handle:
        for raw in handle:
            chrom, start, end, *_, stain = raw.rstrip("\n").split("\t")
            if stain == "acen":
                lo, hi = cen.get(chrom, (10**12, 0))
                cen[chrom] = (min(lo, int(start)), max(hi, int(end)))
    return cen


def somatic_junctions(tumor, normal, windows, p, cytoband=None):
    """Junctions: the junctions of the regions, as dicts with a, b, support."""
    cen = centromeres(cytoband) if cytoband else {}
    tol = p["germline_tol"]
    normal_ends = defaultdict(list)
    for a, b, _ in normal:
        for x, y in ((a, b), (b, a)):
            normal_ends[(x[0], x[1] // tol)].append((x, y))

    def germline(a, b):
        for k in (-1, 0, 1):
            for x, y in normal_ends[(a[0], a[1] // tol + k)]:
                if (abs(x[1] - a[1]) <= tol and x[2] == a[2] and y[0] == b[0]
                        and abs(y[1] - b[1]) <= tol and y[2] == b[2]):
                    return True
        return False

    def usable(end):
        if end[0] not in cen:
            return True
        lo, hi = cen[end[0]]
        return not lo - p["centromere_margin"] <= end[1] <= hi + p["centromere_margin"]

    def small_deletion(a, b):
        return (a[0] == b[0] and a[2] == "R" and b[2] == "L"
                and b[1] - a[1] <= p["max_small_deletion"])

    def inside(end):
        # an R end is the right end of a segment, so it needs sequence of the window to its
        # left (pos > lo); an L end needs sequence to its right (pos < hi)
        c, pos, side = end
        return any(c == w and (lo < pos <= hi if side == "R" else lo <= pos < hi)
                   for w, lo, hi in windows)

    by_sides = defaultdict(list)
    for a, b, n in tumor:
        by_sides[(a[0], a[2], b[0], b[2])].append((a[1], b[1], n))

    def support(a, b):
        best, w = 0, p["cluster_window"]
        for x, y in ((a, b), (b, a)):
            for pa, pb, n in by_sides.get((x[0], x[2], y[0], y[2]), []):
                if abs(pa - x[1]) <= w and abs(pb - y[1]) <= w:
                    best = max(best, n)
        return best

    out = []
    for a, b, n in tumor:
        if n < p["min_support"] or not (usable(a) and usable(b)) or small_deletion(a, b):
            continue
        if germline(a, b) or not (inside(a) and inside(b)):
            continue
        out.append(dict(a=a, b=b, support=support(a, b)))
    out.sort(key=lambda j: (chrom_key(j["a"][0]), j["a"][1]))
    for k, j in enumerate(out):
        j["name"] = f"j{k}"
    return out


# ---------------------------------------------------------------- segments

def cut_segments(windows, junctions):
    """Segments: the segments (dicts name, chrom, cen, tel, window) and the junction ends."""
    segs, seg_at = [], {}
    for i, (c, lo, hi) in enumerate(windows):
        cuts = sorted({lo, hi} | {e[1] for j in junctions for e in (j["a"], j["b"])
                                  if e[0] == c and lo <= e[1] <= hi})
        for a, b in zip(cuts, cuts[1:]):
            segs.append(dict(chrom=c, cen=a, tel=b, window=i))
    for s, name in zip(segs, segment_ids(len(segs))):
        s["name"] = name
        seg_at[(s["chrom"], s["cen"], "L")] = (name, "cen")
        seg_at[(s["chrom"], s["tel"], "R")] = (name, "tel")
    for j in junctions:
        j["end_a"], j["end_b"] = seg_at[j["a"]], seg_at[j["b"]]
    return segs


def chains(segs):
    out = defaultdict(list)
    for s in segs:
        out[s["window"]].append(s["name"])
    return [out[w] for w in sorted(out)]


# ---------------------------------------------------------------- depth

def _bedcov(bam, bins, mapq):
    text = "".join(f"{c}\t{s}\t{e}\n" for c, s, e in bins)
    with tempfile.NamedTemporaryFile("w", suffix=".bed") as handle:
        handle.write(text)
        handle.flush()
        out = pysam.bedcov("-Q", str(mapq), handle.name, bam)
    d = {}
    for ln in out.splitlines():
        c, s, e, n = ln.split("\t")[:4]
        d[(c, int(s))] = int(n) / (int(e) - int(s))
    return d


def _base_depth(bam, c, lo, hi, mapq):
    out = pysam.depth("-a", "-Q", str(mapq), "-r", f"{c}:{lo + 1}-{hi}", bam)
    v = [int(ln.split("\t")[2]) for ln in out.splitlines()]
    return statistics.mean(v) if v else 0.0


def segment_depth(windows, segs, tumor_bam, normal_bam, p):
    """Depth: {segment: corrected tumor depth}."""
    size = p["bin"]
    bins = [(c, s, s + size) for c, lo, hi in windows for s in range(lo, hi - size + 1, size)]
    T, N = _bedcov(tumor_bam, bins, p["depth_mapq"]), _bedcov(normal_bam, bins, p["depth_mapq"])
    med = statistics.median(N.values()) if N else 0
    if med <= 0:
        raise ValueError(f"{normal_bam}: no normal read depth over the regions (MAPQ >= {p['depth_mapq']})")
    good = {k for k, v in N.items() if abs(v - med) <= p["mapp_tol"] * med}
    out = {}
    for s in segs:
        c, lo, hi = s["chrom"], s["cen"], s["tel"]
        trim = min(p["trim"], (hi - lo) // 4)
        ps = [x for (cc, x) in good if cc == c and x >= lo + trim and x + size <= hi - trim]
        if len(ps) >= 3:
            out[s["name"]] = statistics.median(T[(c, x)] / N[(c, x)] * med for x in ps)
        else:
            td = _base_depth(tumor_bam, c, lo, hi, p["depth_mapq"])
            nd = _base_depth(normal_bam, c, lo, hi, p["depth_mapq"])
            out[s["name"]] = td / nd * med if nd > 0 and abs(nd - med) <= 0.25 * med else td
    return out


# ---------------------------------------------------------------- subpaths

class Locus:
    """Segments, chains and edges of the regions, for walking reads."""

    def __init__(self, segs, junctions):
        self.segs = segs
        self.seg = {s["name"]: s for s in segs}
        self.chains = chains(segs)
        self.edges = {}
        for ch in self.chains:
            for a, b in zip(ch, ch[1:]):
                self.edges[frozenset(((a, "tel"), (b, "cen")))] = f"ref:{a}{b}"
        for j in junctions:
            self.edges[frozenset((j["end_a"], j["end_b"]))] = j["name"]
        self.by_chrom = defaultdict(list)
        for s in segs:
            self.by_chrom[s["chrom"]].append(s)

    def parse(self, t):
        return (t[1:], False) if t.startswith("-") else (t, True)

    def step(self, a, b):
        sa, fa = self.parse(a)
        sb, fb = self.parse(b)
        return self.edges.get(frozenset(((sa, "tel" if fa else "cen"),
                                         (sb, "cen" if fb else "tel"))))


def _rc(tokens):
    return [t[1:] if t.startswith("-") else "-" + t for t in reversed(tokens)]


def _canonical(tokens):
    """The reading of a constraint that is written: more forward segments, then the smaller
    as integers (1 2 3 rather than -3 -2 -1)."""
    def key(ts):
        return (-sum(not t.startswith("-") for t in ts), [int(t) for t in ts])
    return min(list(tokens), _rc(tokens), key=key)


def _blocks(read, locus, p):
    recs = [(read.reference_name, read.reference_start + 1, "-" if read.is_reverse else "+",
             read.cigarstring, read.mapping_quality)]
    if read.has_tag("SA"):
        for r in read.get_tag("SA").strip(";").split(";"):
            c, pos, st, cg, mq, _ = r.split(",")
            recs.append((c, int(pos), st, cg, int(mq)))
    out = []
    for c, pos, st, cg, mq in recs:
        ops = [(int(n), o) for n, o in _CIGAR.findall(cg)]
        lc = ops[0][0] if ops[0][1] in "SH" else 0
        rcl = ops[-1][0] if ops[-1][1] in "SH" else 0
        qa = sum(n for n, o in ops if o in "MI=X")
        rl = sum(n for n, o in ops if o in "MDN=X")
        qs = lc if st == "+" else rcl
        # a short alignment that overlaps no segment counts as unaligned read (max_gap in
        # subpaths), so it neither cuts the string nor keeps the blocks around it apart
        off = not any(min(pos + rl, s["tel"]) > max(pos, s["cen"])
                      for s in locus.by_chrom.get(c, []))
        if mq >= p["read_mapq"] and not (off and qa <= p["max_gap"]):
            out.append((qs, qs + qa, c, pos, pos + rl, st))
    out.sort()
    merged = []
    for b in out:
        if merged:
            a = merged[-1]
            if a[2] == b[2] and a[5] == b[5]:
                gap = (b[3] - a[4]) if a[5] == "+" else (a[3] - b[4])
                if -p["merge"] <= gap <= p["merge"]:
                    merged[-1] = (a[0], b[1], a[2], min(a[3], b[3]), max(a[4], b[4]), a[5])
                    continue
        merged.append(b)
    return merged


def _tokens(locus, block):
    _, _, c, rs, re_, st = block
    hits = [s for s in locus.by_chrom.get(c, [])
            if min(re_, s["tel"]) - max(rs, s["cen"]) >= min(300, (s["tel"] - s["cen"]) // 2 + 1)]
    if not hits:
        return None
    hits.sort(key=lambda s: s["cen"], reverse=st == "-")
    return [s["name"] if st == "+" else "-" + s["name"] for s in hits]


def subpaths(locus, windows, bam, p):
    """Subpath constraints: Counter {tokens (written orientation): reads} of one read set."""
    words, seen = Counter(), set()
    with pysam.AlignmentFile(bam) as handle:
        for c, lo, hi in windows:
            for read in handle.fetch(region=f"{c}:{max(1, lo - p['pad'])}-{hi + p['pad']}"):
                # a supplementary record stands for the whole read too (its SA tag lists the
                # other alignments), so a read is kept when any of its alignments is fetched
                if read.flag & 0x104 or read.query_name in seen:
                    continue
                seen.add(read.query_name)
                pieces, cur, prev = [], [], None
                for b in _blocks(read, locus, p):
                    toks = _tokens(locus, b)
                    if toks is None or (prev is not None and b[0] - prev[1] > p["max_gap"]):
                        if cur:
                            pieces.append(cur)
                        cur = []
                    for tok in toks or []:
                        if cur and locus.step(cur[-1], tok) is None:
                            pieces.append(cur)
                            cur = []
                        cur = cur + [tok]
                    prev = b
                if cur:
                    pieces.append(cur)
                for piece in pieces:
                    words[tuple(_canonical(piece))] += 1
    return words


# ---------------------------------------------------------------- fit

def _cn_scale(locus, junctions, depth, cn_files, line_unit, p):
    """Depth per copy on the scale of the copy-number calls (the prior's center)."""
    if not cn_files:
        return line_unit
    hp = []
    for path in cn_files:
        rows = defaultdict(list)
        with open(path, encoding="utf-8") as handle:
            for ln in handle:
                if not ln.strip() or ln.startswith(("#", "track", "browser")):
                    continue
                c, s, e, _, cn = ln.split("\t")[:5]
                rows[c].append((int(s), int(e), float(cn)))
        hp.append(rows)
    bps = {(locus.seg[n]["chrom"], locus.seg[n][end])
           for j in junctions for n, end in (j["end_a"], j["end_b"])}
    cal = []
    for s in locus.segs:
        c, lo, hi = s["chrom"], s["cen"], s["tel"]
        if hi - lo < 50000:
            continue
        total, ok = 0.0, True
        for h in hp:
            cov = 0
            for a, b, cn in h[c]:
                if a < hi and b > lo:
                    ov = min(hi, b) - max(lo, a)
                    cov += ov
                    total += ov * cn / (hi - lo)
            ok &= cov >= 0.9 * (hi - lo)
        cuts = sorted({x for h in hp for a, b, _ in h[c] for x in (a, b)})
        for a, b in zip(cuts, cuts[1:]):
            if a < hi and b > lo and any(pc == c and a + p["tol"] < q < b - p["tol"] for pc, q in bps):
                ok = False
        if ok and total > 0:
            cal.append((min(hi - lo, 300000), depth[s["name"]], total))
    if not cal:
        return line_unit
    return sum(w * d * c for w, d, c in cal) / sum(w * c * c for w, d, c in cal)


def fit(locus, junctions, depth, line_unit, scale, bounds, p):
    """Fit: the best graph (cn, mult, ref, ends) over depths per copy."""
    edge_lb, seg_lb = bounds
    use_bounds = bounds != (Counter(), Counter())
    CH = locus.chains
    JN = [j["name"] for j in junctions]
    JEND = [[j["end_a"], j["end_b"]] for j in junctions]
    support = {j["name"]: j["support"] for j in junctions}
    length = {s["name"]: s["tel"] - s["cen"] for s in locus.segs}
    order = [s["name"] for s in locus.segs]
    ref_name = {(a, "tel"): f"ref:{a}{b}" for ch in CH for a, b in zip(ch, ch[1:])}

    def propagate(mult, heads):
        jd = Counter()
        for (a, b), k in zip(JEND, mult):
            jd[a] += k
            jd[b] += k
        cn, ref, ends = {}, {}, {}
        for ch, c0 in zip(CH, heads):
            first = ch[0]
            e = c0 - jd[(first, "cen")]
            if e < 0:
                return None
            if e:
                ends[f"{first}.cen"] = e
            cn[first] = c0
            for a, b in zip(ch, ch[1:]):
                r = cn[a] - jd[(a, "tel")]
                if r < 0 or (use_bounds and r < edge_lb[ref_name[(a, "tel")]]):
                    return None
                ref[ref_name[(a, "tel")]] = r
                cn[b] = r + jd[(b, "cen")]
            e = cn[ch[-1]] - jd[(ch[-1], "tel")]
            if e < 0:
                return None
            if e:
                ends[f"{ch[-1]}.tel"] = e
        if use_bounds and any(cn[s] < k for s, k in seg_lb.items()):
            return None
        total = sum(ends.values())
        if total < 2 or total % 2:
            return None
        return cn, ref, ends

    def deviance(mult):
        js = [j for j, k in zip(JN, mult) if support[j]]
        m = dict(zip(JN, mult))
        lam = sum(support[j] for j in js) / max(1e-9, sum(max(m[j], 0.25) for j in js))
        return sum(2 * (support[j] * math.log(support[j] / (lam * max(m[j], 0.25)))
                        - (support[j] - lam * max(m[j], 0.25))) for j in js)

    rng = random.Random(p["seed"])
    nj = len(JN)
    lb = [edge_lb[j] if use_bounds else 0 for j in JN] + [0] * len(CH)
    lb = [max(b, 1) if i < nj else b for i, b in enumerate(lb)]
    graphs = {}
    warm = None
    shape = {}  # x -> (cn, ref, ends, deviance) or None; independent of the depth per copy

    def shape_of(x):
        if x not in shape:
            if len(shape) >= 50_000:  # bounded: large loci visit millions of candidates
                shape.clear()
            r = None
            if all(v >= b for v, b in zip(x, lb)):
                res = propagate(x[:nj], x[nj:])
                if res is not None:
                    r = res + (deviance(x[:nj]),)
            shape[x] = r
        return shape[x]

    n = int(round((p["hi"] - p["lo"]) / p["step"])) + 1
    for i in range(n):
        u = line_unit * (p["lo"] + i * p["step"])
        prior = ((u - scale) / (p["prior_sd"] * scale)) ** 2 if p["prior_sd"] > 0 else 0.0
        obs = {s: depth[s] / u for s in order}
        cache = {}

        def J(x):
            x = tuple(x)
            if x not in cache:
                r = None
                sh = shape_of(x)
                if sh is not None:
                    cn, ref, ends, dev = sh
                    d = sum(min(length[s], 300000) / 1e5 * (cn[s] - obs[s]) ** 2 for s in obs)
                    r = (d + dev, {"cn": cn, "ref": ref, "ends": ends})
                cache[x] = r
            return cache[x]

        def descend(x):
            cur = J(x)
            kmax = 3 if len(x) <= 30 else 2
            while True:
                found = None
                for k in range(1, kmax + 1):
                    for idx in itertools.combinations(range(len(x)), k):
                        for sg in itertools.product((-1, 1), repeat=k):
                            y = list(x)
                            for a, s in zip(idx, sg):
                                y[a] += s
                            r = J(y)
                            if r is not None and (cur is None or r[0] < cur[0] - 1e-9) and \
                                    (found is None or r[0] < found[1][0]):
                                found = (y, r)
                    if found:
                        break
                if not found:
                    return x, cur
                x, cur = found

        x0 = list(lb[:nj]) + [max(0, round(obs[ch[0]])) for ch in CH]
        # descend() cannot leave an infeasible start: raise the chain heads until x0 is
        # balanced (more head copies only add intact copies, so this ends)
        while J(x0) is None and max(x0[nj:], default=0) < 1000:
            x0 = x0[:nj] + [h + 1 for h in x0[nj:]]
        starts = [x0] + ([warm] if warm else [])
        for k in range(p["starts"]):
            base = starts[k % len(starts)]
            starts.append([max(b, v + rng.randint(-2, 2)) for v, b in zip(base, lb)])
        best_x, best = None, None
        for s in starts:
            x, r = descend(s)
            if r is not None and (best is None or r[0] < best[0] - 1e-9):
                best_x, best = x, r
        warm = best_x
        for x, r in cache.items():
            if r is None:
                continue
            j = r[0] + prior
            if x not in graphs or j < graphs[x]["J"] - 1e-9:
                graphs[x] = dict(r[1], J=j, unit=u, mult=dict(zip(JN, x[:nj])))
        if len(graphs) > 20000:
            keep = sorted(graphs, key=lambda k: graphs[k]["J"])[:5000]
            graphs = {k: graphs[k] for k in keep}
    if not graphs:
        return None
    return min(graphs.values(), key=lambda g: g["J"])


# ---------------------------------------------------------------- driver

def foldback_flags(junctions, p):
    """{name: 0/1}: a junction is a fold-back when it has the fold-back shape (is_foldback) and is
    not half of a reciprocal inversion: no other fold-back-shaped junction of the opposite side has
    both ends within reciprocal_tol of its ends."""
    def shape(j):
        return is_foldback(j["a"], j["b"], p["foldback_span"])

    def span(j):
        return sorted((j["a"][1], j["b"][1]))

    out = {}
    for j in junctions:
        out[j["name"]] = int(shape(j) and not any(
            o is not j and shape(o) and o["a"][0] == j["a"][0] and o["a"][2] != j["a"][2]
            and all(abs(x - y) <= p["reciprocal_tol"] for x, y in zip(span(o), span(j)))
            for o in junctions))
    return out


def _bounds(locus, words):
    """Lower bounds (edges, segments): every constraint is present in the graph."""
    edge_lb, seg_lb = Counter(), Counter()
    for toks in words:
        used = [locus.step(a, b) for a, b in zip(toks, toks[1:])]
        for e, k in Counter(used).items():
            edge_lb[e] = max(edge_lb[e], k)
        for s, k in Counter(locus.parse(t)[0] for t in toks).items():
            seg_lb[s] = max(seg_lb[s], k)
    return edge_lb, seg_lb


def build_graph(regions, junctions_tsv, normal_junctions_tsv, tumor_bam, normal_bam,
                haploid_depth, cn_files=(), reads=(), use_subpaths=False, cytoband=None,
                header=None, **params):
    """Junctions, segments, depth, subpath constraints and fit. Returns (GraphFile, subpath rows [(string, reads, read_set)])."""
    p = dict(DEFAULTS, **params)
    windows = sorted((parse_region(r) for r in regions), key=lambda w: (chrom_key(w[0]), w[1]))
    junctions = somatic_junctions(read_calls(junctions_tsv), read_calls(normal_junctions_tsv),
                                  windows, p, cytoband)
    segs = cut_segments(windows, junctions)
    locus = Locus(segs, junctions)
    depth = segment_depth(windows, segs, tumor_bam, normal_bam, p)

    rows, all_words = [], set()
    for name, bam in reads:
        words = subpaths(locus, windows, bam, p)
        all_words |= set(words)
        rows += [(" ".join(w), n, name) for w, n in words.most_common()]

    scale = _cn_scale(locus, junctions, depth, list(cn_files), haploid_depth, p)
    bounds = _bounds(locus, all_words) if use_subpaths else (Counter(), Counter())
    best = fit(locus, junctions, depth, haploid_depth, scale, bounds, p)
    if best is None:
        raise ValueError("no balanced graph satisfies the subpath constraints")
    cn, mult, ref, free = best["cn"], best["mult"], best["ref"], best["ends"]

    side = {"cen": "L", "tel": "R"}

    def end(name, e):
        return f"{name}.{side[e]}"

    flags = foldback_flags(junctions, p)

    graph = GraphFile(header=dict(header or {}, haploid_depth=f"{best['unit']:.4f}"))
    graph.segments = [GraphSegment(s["name"], s["chrom"], s["cen"], s["tel"], "+", cn[s["name"]])
                      for s in segs]
    for ch in locus.chains:
        for x, y in zip(ch, ch[1:]):
            if ref[f"ref:{x}{y}"] > 0:
                graph.reference.append((end(x, "tel"), end(y, "cen"), ref[f"ref:{x}{y}"]))
    for j in junctions:
        if mult[j["name"]] > 0:
            graph.junctions.append((j["name"], end(*j["end_a"]), end(*j["end_b"]),
                                    mult[j["name"]], flags[j["name"]]))
    order = {s["name"]: k for k, s in enumerate(segs)}
    for e, k in sorted(free.items(), key=lambda x: (order[x[0].split(".")[0]], x[0])):
        name, which = e.split(".")
        graph.free_ends.append((end(name, which), k))
    graph.validate()
    return graph, rows
