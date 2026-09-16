"""The ``O(2^m n^4)`` algorithm for the nested BFB distance.

Every node of the NBFB graph is represented by O(1) indices into the
target ``t`` (paper, appendix "An O(2^m n^4)-Time Algorithm"):

* an empty-tail node ``E(phi)`` is the prefix ``t[1..phi]``;
* a window node ``W(phi, a, b)`` is ``t[1..phi]`` followed by the reverse
  complement of ``t[a..b]``, with matched prefix exactly ``phi``;
* a source node ``T(phi0, x)`` is ``t[1..phi0]`` followed by source ``x``,
  produced by a t-BFB cycle, with ``1 <= j < |x|`` symbols of ``x`` matched.

Successors are computed from longest-common-extension tables and an
occurrence table, and the successors of window and source nodes are
processed as intervals of tuples over an interval union-find structure, so
every tuple is visited once in amortized O(1) time. All indices in this
module are 1-based to mirror the paper; ``Sequence`` slicing converts at the
boundary.

The search records a parent pointer per tuple, so it returns one shortest
history rather than every shortest history. Use ``NBFB_dist`` for the
complete minimum-path graph.
"""

from collections import deque
import math

from bfb_dist.bfb_graph import BFBGraph
from bfb_dist.history import (
    Operation,
    State,
    start_state,
    validate_sources,
)


def _lce_table(y, z):
    """Return the 1-based table ``LCE[i][j]`` of ``y[i..]`` versus ``z[j..]``.

    Rows ``0`` and ``len(y) + 1`` and columns ``0`` and ``len(z) + 1`` are
    zero, so callers may query one position past either end.
    """
    ly, lz = len(y), len(z)
    table = [[0] * (lz + 2) for _ in range(ly + 2)]
    for i in range(ly, 0, -1):
        yi = y[i - 1]
        row, next_row = table[i], table[i + 1]
        for j in range(lz, 0, -1):
            if yi == z[j - 1]:
                row[j] = 1 + next_row[j + 1]
    return table


class _NextUnvisited:
    """Interval union-find over ``lo..hi`` supporting next-unvisited queries.

    ``next(v)`` returns the smallest unvisited value at or above ``v``, or
    ``hi + 1`` when there is none. ``remove(v)`` marks ``v`` visited. Both
    take amortized near-constant time through path compression.
    """

    def __init__(self, lo, hi):
        self.lo = lo
        self.hi = hi
        # parent[v - lo] == v while v is unvisited; hi + 1 is a sentinel.
        self.parent = list(range(lo, hi + 2))

    def next(self, v):
        if v < self.lo:
            v = self.lo
        if v > self.hi:
            return self.hi + 1
        lo, parent = self.lo, self.parent
        root = v
        while parent[root - lo] != root:
            root = parent[root - lo]
        while parent[v - lo] != root:
            parent[v - lo], v = root, parent[v - lo]
        return root

    def remove(self, v):
        self.parent[v - self.lo] = v + 1

    def is_unvisited(self, v):
        return self.parent[v - self.lo] == v


class _FastSearch:
    """Breadth-first search over index tuples for one instance."""

    def __init__(self, sources, target):
        self.sources = sources
        self.target = target
        self.n = n = len(target)
        self.m = len(sources)
        t = target.segments
        tbar = target.reverse_complement().segments
        self.lce_tt = _lce_table(t, t)
        self.lce_bt = _lce_table(tbar, t)
        self.lce_tb = _lce_table(t, tbar)
        self.source_len = [len(x) for x in sources]
        self.lce_xt = [_lce_table(x.segments, t) for x in sources]
        self.lce_tx = [_lce_table(t, x.segments) for x in sources]
        self.lce_txbar = [
            _lce_table(t, x.reverse_complement().segments) for x in sources
        ]
        self.z_len = {"bar": n}
        self.lce_tz = {"bar": self.lce_tb}
        for index in range(self.m):
            self.z_len[index] = self.source_len[index]
            self.lce_tz[index] = self.lce_tx[index]
        self._occ = {}
        self._minb = {}

        # Tables for t-BFB cycles that keep the matched prefix. Such a cycle
        # appends a source y to a nonempty tail u and is admitted when the
        # produced tail u y is foldable, i.e. when the reverse complement
        # ybar ubar occurs in t. ``xbar_pos[y]`` lists the start positions of
        # ybar in t; ``ext_w[y][a]`` is the longest match of t[a..] after such
        # an occurrence, with its position; ``ext_t[y][x]`` is the longest
        # match of xbar after such an occurrence, with its position.
        self.xbar_pos = []
        for index in range(self.m):
            length = self.source_len[index]
            table = self.lce_txbar[index]
            self.xbar_pos.append(
                [p for p in range(1, n - length + 2) if table[p][1] >= length]
            )
        self.ext_w = []
        for y in range(self.m):
            ly = self.source_len[y]
            row = [(0, 0)] * (n + 2)
            for a in range(1, n + 1):
                best = (0, 0)
                for p in self.xbar_pos[y]:
                    value = self.lce_tt[p + ly][a]
                    if value > best[0]:
                        best = (value, p)
                row[a] = best
            self.ext_w.append(row)
        self.ext_t = [[(0, 0)] * self.m for _ in range(self.m)]
        for y in range(self.m):
            ly = self.source_len[y]
            for x in range(self.m):
                best = (0, 0)
                for p in self.xbar_pos[y]:
                    value = self.lce_txbar[x][p + ly][1]
                    if value > best[0]:
                        best = (value, p)
                self.ext_t[y][x] = best

        self.dist = {}
        self.parent = {}
        self.e_family = {}
        self.w_family = {}
        self.queue = deque()

    # -- tables -----------------------------------------------------------

    def occ(self, zkey, i, lam, s):
        """Return ``(OCC_z(i, lam, s), POS_z(i, lam, s))``."""
        key = (zkey, i, lam, s)
        cached = self._occ.get(key)
        if cached is not None:
            return cached
        lce_tt, lce_tz = self.lce_tt, self.lce_tz[zkey]
        best, pos = -1, i
        for p in range(1, self.n - lam + 2):
            if lce_tt[p][i] >= lam:
                value = lce_tz[p + lam][s]
                if value > best:
                    best, pos = value, p
        result = (best, pos)
        self._occ[key] = result
        return result

    def minb(self, zkey, phi, s, mu):
        """Return ``(MINB_z(phi, s, mu), POS at MINB)``; MINB is ``phi`` if none."""
        key = (zkey, phi, s, mu)
        cached = self._minb.get(key)
        if cached is not None:
            return cached
        # OCC_z(beta + 1, phi - beta, s) is nondecreasing in beta.
        lo, hi = 0, phi - 1
        answer, pos = phi, None
        while lo <= hi:
            mid = (lo + hi) // 2
            value, p = self.occ(zkey, mid + 1, phi - mid, s)
            if value >= mu:
                answer, pos = mid, p
                hi = mid - 1
            else:
                lo = mid + 1
        result = (answer, pos)
        self._minb[key] = result
        return result

    # -- visiting ---------------------------------------------------------

    def _e_family(self, mask):
        family = self.e_family.get(mask)
        if family is None:
            family = self.e_family[mask] = _NextUnvisited(0, self.n)
        return family

    def _w_family(self, mask, phi, b):
        key = (mask, phi, b)
        family = self.w_family.get(key)
        if family is None:
            family = self.w_family[key] = _NextUnvisited(1, self.n)
        return family

    def visit(self, mask, node, parent, operation):
        """Record ``node`` at distance ``dist[parent] + 1`` if new."""
        kind = node[0]
        if kind == "E":
            family = self._e_family(mask)
            if not family.is_unvisited(node[1]):
                return
            family.remove(node[1])
        elif kind == "W":
            family = self._w_family(mask, node[1], node[3])
            if not family.is_unvisited(node[2]):
                return
            family.remove(node[2])
        else:
            if (mask, node) in self.dist:
                return
        key = (mask, node)
        self.dist[key] = self.dist[parent] + 1
        self.parent[key] = (parent, operation)
        self.queue.append(key)

    def _visit_interval(self, family, make_node, lo, hi, mask, parent, make_op):
        """Visit every unvisited value of ``family`` in ``lo..hi``."""
        v = family.next(lo)
        while v <= hi:
            self.visit(mask, make_node(v), parent, make_op(v))
            v = family.next(v + 1)

    # -- expansion --------------------------------------------------------

    def expand_empty(self, mask, phi, key):
        n = self.n
        lce_bt = self.lce_bt

        for index in range(self.m):
            if mask >> index & 1:
                continue
            length = self.source_len[index]
            j = min(length, self.lce_xt[index][1][phi + 1])
            if j >= 1:
                node = ("E", phi + length) if j == length else ("T", phi, index)
            elif self.xbar_pos[index]:
                # The prefix is kept; the tail x is foldable as xbar occurs.
                p = self.xbar_pos[index][0]
                node = ("W", phi, p, p + length - 1)
            else:
                continue
            self.visit(
                mask | 1 << index, node, key, Operation.translocation(index)
            )

        for beta in range(phi):
            # p-BFB: q = t[beta+1..phi].
            M = min(phi - beta, lce_bt[n - phi + 1][phi + 1])
            if M >= 1:
                node = (
                    ("E", phi + M)
                    if M == phi - beta
                    else ("W", phi + M, beta + 1, phi - M)
                )
                self.visit(mask, node, key, Operation("p", beta, phi))

            for gamma in range(beta + 1, phi):
                # i-BFB: q = t[beta+1..gamma], r = t[gamma+1..phi].
                M = min(gamma - beta, lce_bt[n - gamma + 1][phi + 1])
                if M >= 1:
                    node = (
                        ("E", phi + M)
                        if M == gamma - beta
                        else ("W", phi + M, beta + 1, gamma - M)
                    )
                    self.visit(mask, node, key, Operation("i", beta, gamma))

                # f-BFB: q = t[beta+1..c], r = t[c+1..phi], c = gamma.
                c = gamma
                M = min(phi - beta, lce_bt[n - phi + 1][c + 1])
                if c + M > phi:
                    node = (
                        ("E", c + M)
                        if M == phi - beta
                        else ("W", c + M, beta + 1, phi - M)
                    )
                    self.visit(mask, node, key, Operation("f", beta, c))

    def expand_tail(self, mask, node, key):
        n = self.n
        lce_bt = self.lce_bt
        if node[0] == "W":
            _, phi, a, b = node
            lam = b - a + 1
            lce_row_for_c = lambda c: self.lce_tt[c + 1][a]
            zkey, s = "bar", n - b + 1
        else:
            _, phi0, index = node
            length = self.source_len[index]
            j = min(length, self.lce_xt[index][1][phi0 + 1])
            phi = phi0 + j
            lam = length - j
            lce_txbar = self.lce_txbar[index]
            lce_row_for_c = lambda c: lce_txbar[c + 1][1]
            zkey, s = index, j + 1

        for c in range(1, phi + 1):
            M = min(lam, lce_row_for_c(c))
            if M == 0:
                continue

            if M == lam:
                L = lce_bt[n - phi + 1][c + lam + 1]
                # Full absorption, empty-tail successors E(phi') with
                # beta = phi - (phi' - c - lam).
                lo, hi = lam + phi + 1, c + lam + min(phi, L)
                if lo <= hi:
                    self._visit_interval(
                        self._e_family(mask),
                        lambda v: ("E", v),
                        lo,
                        hi,
                        mask,
                        key,
                        lambda v, c=c: Operation("f", phi - (v - c - lam), c),
                    )
                # Full absorption, window successors W(c+lam+L, a', phi-L)
                # with beta = a' - 1.
                if c + lam + L > phi:
                    hi_a = min(c, phi - L)
                    if hi_a >= 1:
                        new_phi, new_b = c + lam + L, phi - L
                        self._visit_interval(
                            self._w_family(mask, new_phi, new_b),
                            lambda v: ("W", new_phi, v, new_b),
                            1,
                            hi_a,
                            mask,
                            key,
                            lambda v, c=c: Operation("f", v - 1, c),
                        )
                continue

            # Partial absorption.
            if c + M <= phi:
                continue
            mu = lam - M
            beta_min, P = self.minb(zkey, phi, s, mu)
            if beta_min > c - 1:
                continue
            new_phi = c + M
            new_b = P - beta_min + phi + mu - 1
            self._visit_interval(
                self._w_family(mask, new_phi, new_b),
                lambda v: ("W", new_phi, v, new_b),
                P,
                P + c - 1 - beta_min,
                mask,
                key,
                lambda v, c=c: Operation("f", beta_min + (v - P), c),
            )

        # t-BFB cycles that keep the matched prefix: appending y to the tail
        # u is admitted when ybar ubar occurs in t, and the successor is the
        # window node whose window spells ybar ubar.
        for y in range(self.m):
            if mask >> y & 1:
                continue
            ly = self.source_len[y]
            if node[0] == "W":
                L, p = self.ext_w[y][a]
            else:
                L, p = self.ext_t[y][index]
            if L >= lam:
                self.visit(
                    mask | 1 << y,
                    ("W", phi, p, p + ly + lam - 1),
                    key,
                    Operation.translocation(y),
                )

    # -- search -----------------------------------------------------------

    def run(self):
        """Return the goal key if reached, else ``None``."""
        start = (0, ("E", 0))
        self.dist[start] = 0
        self.parent[start] = None
        self._e_family(0).remove(0)
        self.queue.append(start)

        while self.queue:
            key = self.queue.popleft()
            mask, node = key
            if node == ("E", self.n):
                return key
            if node[0] == "E":
                self.expand_empty(mask, node[1], key)
            else:
                self.expand_tail(mask, node, key)
        return None

    def to_state(self, key):
        """Return the :class:`State` a tuple represents."""
        mask, node = key
        t = self.target
        if node[0] == "E":
            sequence = t[:node[1]]
        elif node[0] == "W":
            _, phi, a, b = node
            sequence = t[:phi] + t[a - 1:b].reverse_complement()
        else:
            _, phi0, index = node
            sequence = t[:phi0] + self.sources[index]
        return State(sequence, mask)

    def history(self, goal):
        """Return the list of ``(state, operation)`` pairs ending at ``goal``."""
        steps = []
        key = goal
        while self.parent[key] is not None:
            parent, operation = self.parent[key]
            steps.append((parent, operation))
            key = parent
        steps.reverse()
        return steps


def NBFB_fast_dist(sources, target):
    """Compute the nested BFB distance in ``O(2^m n^4)`` time.

    Returns a :class:`~bfb_dist.bfb_graph.BFBGraph` containing exactly one
    shortest nested history, or only the start state when the target is
    unreachable. The reconstructed history is verified by replaying its
    operations on sequences.

    Raises:
        RuntimeError: If the replayed history does not reach the target,
            which indicates an inconsistency between the index-tuple rules
            and the sequence-level definitions.
    """
    sources = validate_sources(sources)
    start = start_state()
    if len(target) == 0:
        return BFBGraph({start: {}}, sources, target, 0)

    search = _FastSearch(sources, target)
    goal = search.run()
    if goal is None:
        return BFBGraph({start: {}}, sources, target, math.inf)

    adjacency = {}
    state = start
    for key, operation in search.history(goal):
        expected = search.to_state(key)
        if state != expected:
            raise RuntimeError(
                f"index tuple {key} represents {expected} but the replayed "
                f"history reached {state}"
            )
        child = state.apply(operation, sources)
        adjacency[state] = {child: operation}
        state = child
    if state.sequence != target:
        raise RuntimeError(
            f"replayed history reached {state} instead of the target"
        )
    return BFBGraph(adjacency, sources, target, len(adjacency))
