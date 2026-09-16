"""Candidate-source extension experiment (2026-09-16).

Sources are every run of consecutive segments of an input genome, in both
orientations. Search: nested class with the relaxed t-BFB rule, sources used
at most once, two consecutive t-BFBs may not append genome-adjacent blocks
(the merged block dominates), iterative deepening on the number of sources
used. Prints the (sources used, distance) frontier rows. Run under
`ulimit -v 4000000`; COLO829 with inserts takes ~5 min for the second row.
"""
from collections import deque

def rc(s): return s[::-1].swapcase()

def blocks_of(genomes):
    """Every run of consecutive segments of an input genome, both orientations.
    A block is (genome id, start, end, orient); its string is derived."""
    out = []
    for gi, g in enumerate(genomes):
        for i in range(len(g)):
            for j in range(i + 1, len(g) + 1):
                out.append((gi, i, j, +1, g[i:j]))
                out.append((gi, i, j, -1, rc(g[i:j])))
    return out

def adjacent(y, x):
    """x continues y in the same genome and orientation (so y+x is a block)."""
    if y is None or y[0] != x[0] or y[3] != x[3]: return False
    return (y[2] == x[1]) if y[3] == +1 else (x[2] == y[1])

def lcp(v, t):
    k = 0
    while k < len(v) and k < len(t) and v[k] == t[k]: k += 1
    return k

def foldable(v, t):
    tail = v[lcp(v, t):]
    return tail == "" or rc(tail) in t

def search(genomes, t, max_sources, max_states=3_000_000):
    B = blocks_of(genomes); n = len(t)
    start = ("", frozenset(), None)
    dist = {start: 0}; parent = {start: None}; q = deque([start]); goals = {}
    while q:
        v, U, last = q.popleft(); d = dist[(v, U, last)]
        if v == t:
            key = U
            if key not in goals: goals[key] = (d, (v, U, last))
            continue
        phi = lcp(v, t)
        succ = []
        # folds
        for x in range(len(v)):
            for y in range(x + 1, len(v) + 1):
                Q, R = v[x:y], v[y:]
                kids = [(v + rc(Q), ("p" if not R else "i", x, y))]
                if R: kids.append((v[:y] + rc(R) + rc(Q), ("f", x, y)))
                for w, op in kids:
                    if len(w) <= n and lcp(w, t) > phi and foldable(w, t):
                        succ.append(((w, U, None), op))
        # t-BFBs
        for bid, b in enumerate(B):
            if bid in U or adjacent(last, b) or len(U) >= max_sources: continue
            w = v + b[4]
            if len(w) > n: continue
            if lcp(w, t) > phi or foldable(w, t):
                succ.append(((w, U | {bid}, b), ("t", b[4])))
        for s, op in succ:
            if s not in dist:
                dist[s] = d + 1; parent[s] = ((v, U, last), op); q.append(s)
                if len(dist) > max_states:
                    import collections
                    by_len = collections.Counter(len(k[0]) for k in dist)
                    by_u = collections.Counter(len(k[1]) for k in dist)
                    sample = [ (k[0], sorted(B[i][4] for i in k[1]), k[2] and k[2][4]) for k in list(dist)[-8:] ]
                    raise MemoryError(f"> {max_states} states; by_len={sorted(by_len.items())}; by_|U|={sorted(by_u.items())}; sample={sample}")
    return B, dist, parent, goals

def ops_of(parent, s):
    ops = []; cur = s
    while parent[cur] is not None:
        prev, op = parent[cur]; ops.append(op); cur = prev
    return ops[::-1]

def report(name, genomes, t, extra_rows=1):
    """Iterative deepening on the number of sources: k = 1, 2, ... until the
    first finite row, then ``extra_rows`` more rows of the frontier."""
    goals = {}; kmax = None
    for k in range(1, len(blocks_of(genomes)) + 1):
        t0 = time.time()
        try:
            B, dist, parent, g = search(genomes, t, k)
        except MemoryError as e:
            print(f"{name}: k={k} aborted ({str(e)[:200]})"); return
        print(f"{name}: k<={k} states={len(dist)} time={time.time()-t0:.1f}s goals={len(g)}")
        for U, (d, st) in g.items():
            if U not in goals: goals[U] = (d, ops_of(parent, st))
        if g and kmax is None: kmax = k + extra_rows
        if kmax is not None and k >= kmax: break
    best = {}
    for U, (d, ops) in goals.items():
        j = len(U)
        if j not in best or d < best[j][0]: best[j] = (d, ops)
    for j in sorted(best):
        d, ops = best[j]
        kinds = "".join(o[0] for o in ops)
        runs = sum(1 for i, k in enumerate(kinds) if k == "t" and (i == 0 or kinds[i-1] != "t"))
        used = [o[1] for o in ops if o[0] == "t"]
        print(f"  sources used={j} distance={d} folds={d-j} t-runs={runs} kinds={kinds} used={used}")
    sys.stdout.flush()

report("HG008 chr19 (germline ZABCDEFGP)", ["ZABCDEFGP"], "ZABCFEFGPgfefcbaB")
report("COLO829 chr3p, no inserts (arm only)", ["ZYCBXAEW"], "ZYCBXAEWeaxbcYCcyCBXAEweaXA")
report("COLO829 chr3p with inserts (arm, chr6 UV, chr10 N, chr12 M)", ["ZYCBXAEW", "UV", "N", "M"], "ZYCBXAEWeaxbcYCUVucyCBXAEweamnXA")
