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
    start_state,
    validate_sources,
)


def _is_foldable(sequence, target):
    """Return whether the tail of ``sequence`` beyond its LCP is foldable."""
    tail = sequence[sequence.lcp_length(target):]
    return len(tail) == 0 or tail.reverse_complement().is_infix_of(target)


def _truncate_G(G, start, ends):
    """Return adjacency rows for nodes with a recorded path to an end state.

    The search explores forward from ``start``, so every key of ``G`` is
    reachable from ``start`` by construction; the missing property is
    reachability *to* a state in ``ends``. This walks backward from the end
    states over the reversed edges to find it.
    """
    parents_of = defaultdict(list)
    for node, children in G.items():
        for child in children:
            parents_of[child].append(node)

    reachable = set(ends)
    pending = deque(ends)
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


def _successors(node, sources, target, foldable, max_translocations):
    """Yield ``(child, operation)`` pairs of one BFB cycle.

    Children of p-, i-, and f-BFB cycles must strictly extend the matched
    prefix and, when ``foldable`` is true, have a foldable tail. Children of
    t-BFB cycles may keep the matched prefix unchanged; a child that extends
    the prefix is exempt from foldability, as the appended source is given
    rather than produced by a fold-back, whereas a child that keeps the
    prefix must have a foldable tail when ``foldable`` is true. Children
    longer than the target are discarded, as no cycle shortens a sequence,
    and so are t-BFB children beyond ``max_translocations`` cycles.
    """
    sequence = node.sequence
    phi = sequence.lcp_length(target)
    n = len(target)

    for child, x, y, is_forward_deletion in sequence.successors():
        if len(child) > n or child.lcp_length(target) <= phi:
            continue
        if foldable and not _is_foldable(child, target):
            continue
        operation = Operation.fold_back(
            x, y, is_forward_deletion, len(sequence)
        )
        yield State(child, node.translocations), operation

    if node.translocations >= max_translocations:
        return
    for index, source in enumerate(sources):
        child = sequence + source
        if len(child) > n:
            continue
        if child.lcp_length(target) <= phi and foldable:
            if not _is_foldable(child, target):
                continue
        yield (
            State(child, node.translocations + 1),
            Operation.translocation(index),
        )


def _bounded_search(sources, target, foldable, max_translocations):
    """Breadth-first search with at most ``max_translocations`` t-BFB cycles.

    Returns the layered adjacency, the start state, the end states of the
    first layer reaching the target, and that layer's index.
    """
    start = start_state()
    pending = deque([start])
    G = defaultdict(dict)
    cost = {start: 0}
    bfb_distance = math.inf
    ends = []

    while pending:
        node = pending.popleft()
        if cost[node] >= bfb_distance:
            break

        new_cost = cost[node] + 1
        for child, operation in _successors(
            node, sources, target, foldable, max_translocations
        ):
            if child not in cost:
                cost[child] = new_cost
                pending.append(child)
                if child.sequence == target:
                    bfb_distance = new_cost
                    ends.append(child)

            if cost[child] == new_cost:
                # Equal sources at different indices yield parallel edges;
                # the first index is kept.
                G[node].setdefault(child, operation)

    return G, start, ends, bfb_distance


def _search(sources, target, foldable):
    """Infer every history in the selected class that is minimum first in
    the number of t-BFB cycles and then in the total number of cycles.

    The search deepens the bound on t-BFB cycles from one upward; the
    first bound admitting a history is the minimum number of t-BFB cycles,
    and every history found under that bound uses exactly that many.
    Every t-BFB cycle appends a nonempty source, so no history reaching the
    target has more t-BFB cycles than the target has segments.
    """
    sources = validate_sources(sources)
    start = start_state()
    if len(target) == 0:
        return BFBGraph({start: {}}, sources, target, 0)

    for max_translocations in range(1, len(target) + 1):
        G, start, ends, bfb_distance = _bounded_search(
            sources, target, foldable, max_translocations
        )
        if ends:
            adjacency = _truncate_G(G, start, ends)
            return BFBGraph(adjacency, sources, target, bfb_distance)
    return BFBGraph({start: {}}, sources, target, math.inf)


def NBFB_dist(sources, target):
    """Infer every minimum-length nested history over ``sources`` to ``target``."""
    return _search(sources, target, foldable=True)


def LCP_monotone_dist(sources, target):
    """Infer every minimum-length LCP-monotone history over ``sources``."""
    return _search(sources, target, foldable=False)
