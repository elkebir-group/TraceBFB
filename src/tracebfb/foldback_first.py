"""BFB distance under the fold-back-first objective (backward search).

A history is ranked lexicographically by the number of fold-back-shaped
t-BFB cycles, then the number of t-BFB cycles, then the total number of
cycles. A t-BFB cycle is fold-back-shaped when its join ``u -> v``, from the
last segment ``u`` of the sequence to the first segment ``v`` of the appended
source, has both segments on one chromosome in opposite orientations, with
breakpoints at most ``span`` bp apart. Such a join is what a fold-back cycle
produces, so the objective prefers histories that make fold-back junctions
with fold-back cycles. ``cpp/nbfb --foldback-first`` is the nested
counterpart.
"""

import heapq
import itertools

from tracebfb.distance import _predecessors
from tracebfb.graph_file import FOLDBACK_SPAN
from tracebfb.history import State, validate_sources
from tracebfb.sequence import Sequence


def foldback_join(u, v, span=FOLDBACK_SPAN):
    """Return whether the join ``u -> v`` is fold-back-shaped."""
    if not (hasattr(u, "chr") and hasattr(v, "chr")):
        return False
    if u.chr != v.chr or u.complement == v.complement:
        return False
    return span is None or abs(u.end - v.start) <= span


def foldback_first_BFB_dist(sources, target, span=FOLDBACK_SPAN):
    """Return one optimal history under the fold-back-first objective.

    The result is ``None`` when no history reaches ``target``, and otherwise
    a dict with the cost components ``foldback_t``, ``t``, and ``total``, the
    cycle kinds ``types`` (``T`` marks a fold-back-shaped t-BFB cycle), and
    the operations ``ops``. ``span=None`` drops the distance test.
    """
    sources = validate_sources(sources)
    bound = len(target) + 1
    goal = Sequence()
    best = {target: (0, 0, 0)}
    parent = {target: None}
    tie = itertools.count()
    heap = [((0, 0, 0), next(tie), target)]
    done = set()
    while heap:
        _, _, sequence = heapq.heappop(heap)
        if sequence in done:
            continue
        done.add(sequence)
        if sequence == goal:
            break
        f, t, n = best[sequence]
        for predecessor, operation in _predecessors(
            State(sequence, 0), sources, bound
        ):
            previous = predecessor.sequence
            shaped = (
                operation.kind == "t"
                and len(previous) > 0
                and foldback_join(previous[-1], sources[operation.source][0], span)
            )
            cost = (f + shaped, t + (operation.kind == "t"), n + 1)
            if previous not in best or cost < best[previous]:
                best[previous] = cost
                parent[previous] = (sequence, operation, shaped)
                heapq.heappush(heap, (cost, next(tie), previous))
    if goal not in best:
        return None

    # The search runs backward, so following parents from the empty
    # sequence lists the operations in forward order.
    operations = []
    current = goal
    while parent[current] is not None:
        current, operation, shaped = parent[current]
        operations.append((operation, shaped))
    f, t, n = best[goal]
    return {
        "foldback_t": f,
        "t": t,
        "total": n,
        "types": "".join(
            ("T" if shaped else "t") if op.kind == "t" else op.kind
            for op, shaped in operations
        ),
        "ops": ";".join(str(op) for op, _ in operations),
    }
