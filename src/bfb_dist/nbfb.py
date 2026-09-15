"""Nested and LCP-monotone deletion-aware BFB distance (forward search).

``NBFB_dist`` is the breadth-first search over the NBFB graph of the paper
(Algorithm NBFB-Distance, ``O(2^m n^6)`` time). ``LCP_monotone_dist`` is the
same search without the foldability test, over the exponentially large
LCP-monotone graph of the appendix.
"""

from collections import defaultdict, deque
import math

from bfb_dist.bfb_graph import BFBGraph
from bfb_dist.history import (
    Operation,
    State,
    end_state,
    start_state,
    total_source_length,
    validate_sources,
)


def _is_foldable(sequence, target):
    """Return whether the tail of ``sequence`` beyond its LCP is foldable."""
    tail = sequence[sequence.lcp_length(target):]
    return len(tail) == 0 or tail.reverse_complement().is_infix_of(target)


def _truncate_G(G, start, end):
    """Return adjacency rows for nodes with a recorded path to ``end``.

    The search explores forward from ``start``, so every key of ``G`` is
    reachable from ``start`` by construction; the missing property is
    reachability *to* ``end``. This walks backward from ``end`` over the
    reversed edges to find it.
    """
    parents_of = defaultdict(list)
    for node, children in G.items():
        for child in children:
            parents_of[child].append(node)

    reachable = {end}
    pending = deque([end])
    while pending:
        node = pending.popleft()
        for parent in parents_of[node]:
            if parent not in reachable:
                reachable.add(parent)
                pending.append(parent)

    adjacency = {
        node: {
            child: operation
            for child, operation in children.items()
            if child in reachable
        }
        for node, children in G.items()
        if node in reachable
    }
    return adjacency or {start: {}}


def _successors(node, sources, target, foldable):
    """Yield ``(child, operation)`` pairs that extend the matched prefix.

    Children of p-, i-, and f-BFB cycles must have a foldable tail when
    ``foldable`` is true. Children of t-BFB cycles are exempt, as the
    appended source is given rather than produced by a fold-back.
    """
    sequence = node.sequence
    phi = sequence.lcp_length(target)

    for child, x, y, is_forward_deletion in sequence.successors():
        if child.lcp_length(target) <= phi:
            continue
        if foldable and not _is_foldable(child, target):
            continue
        operation = Operation.fold_back(
            x, y, is_forward_deletion, len(sequence)
        )
        yield State(child, node.used), operation

    for index, source in enumerate(sources):
        if node.has_used(index):
            continue
        child = sequence + source
        if child.lcp_length(target) <= phi:
            continue
        yield State(child, node.used | 1 << index), Operation.translocation(index)


def _search(sources, target, foldable):
    """Infer every minimum-length history in the selected class."""
    sources = validate_sources(sources)
    start, end = start_state(), end_state(sources, target)
    if start == end:
        return BFBGraph({start: {}}, sources, target, 0)
    if total_source_length(sources) > len(target):
        return BFBGraph({start: {}}, sources, target, math.inf)

    pending = deque([start])
    G = defaultdict(dict)
    cost = {start: 0}
    bfb_distance = math.inf

    while pending:
        node = pending.popleft()
        if cost[node] >= bfb_distance:
            break

        new_cost = cost[node] + 1
        for child, operation in _successors(node, sources, target, foldable):
            if child not in cost:
                cost[child] = new_cost
                pending.append(child)
                if child == end:
                    bfb_distance = new_cost

            if cost[child] == new_cost:
                G[node][child] = operation

    adjacency = _truncate_G(G, start, end)
    return BFBGraph(adjacency, sources, target, bfb_distance)


def NBFB_dist(sources, target):
    """Infer every minimum-length nested history over ``sources`` to ``target``."""
    return _search(sources, target, foldable=True)


def LCP_monotone_dist(sources, target):
    """Infer every minimum-length LCP-monotone history over ``sources``."""
    return _search(sources, target, foldable=False)
