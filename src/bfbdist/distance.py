"""Deletion-aware BFB distance over a set of sources (backward search)."""

from collections import defaultdict, deque
import math

from bfbdist.bfb_graph import BFBGraph
from bfbdist.history import Operation, State, start_state, validate_sources


def _relabel_forward(G, start_backward, total):
    """Return forward adjacency rows reachable from the empty-sequence state.

    The backward search counts the t-BFB cycles it has reversed so far.
    Every history found has ``total`` t-BFB cycles, so a node reached after
    reversing ``c`` of them has ``total - c`` applied in the forward
    direction. The result uses that forward count as the state label.
    """

    def forward(node):
        return State(node.sequence, total - node.translocations)

    adjacency = {}
    seen = {start_backward}
    pending = [start_backward]
    while pending:
        node = pending.pop()
        children = G.get(node)
        if not children:
            continue
        row = adjacency.setdefault(forward(node), {})
        for child, operation in children.items():
            row[forward(child)] = operation
            if child not in seen:
                seen.add(child)
                pending.append(child)
    return adjacency or {start_state(): {}}


def _predecessors(node, sources, max_translocations):
    """Yield ``(predecessor, operation)`` pairs producing ``node``.

    Reversed p-, i-, and f-BFB cycles keep the count of reversed t-BFB
    cycles. A reversed t-BFB cycle strips a source that is a suffix of the
    sequence and raises the count, up to ``max_translocations``.
    """
    sequence = node.sequence
    for predecessor, x, y, is_forward_deletion in sequence.previous_bfb():
        operation = Operation.fold_back(
            x, y, is_forward_deletion, len(predecessor)
        )
        yield State(predecessor, node.translocations), operation

    if node.translocations >= max_translocations:
        return
    for index, source in enumerate(sources):
        predecessor = sequence.strip_suffix(source)
        if predecessor is not None:
            yield (
                State(predecessor, node.translocations + 1),
                Operation.translocation(index),
            )


def _bounded_search(sources, target, max_translocations):
    """Backward breadth-first search reversing at most ``max_translocations``
    t-BFB cycles.

    Returns the layered graph, the empty-sequence state of the first layer
    reaching it (``None`` if none does), and that layer's index.
    """
    end = State(target, 0)
    pending = deque([end])
    G = defaultdict(dict)
    cost = {end: 0}
    bfb_distance = math.inf
    start_backward = None

    while pending:
        node = pending.popleft()
        if cost[node] >= bfb_distance:
            break

        new_cost = cost[node] + 1
        for predecessor, operation in _predecessors(
            node, sources, max_translocations
        ):
            if predecessor not in cost:
                cost[predecessor] = new_cost
                pending.append(predecessor)
                if len(predecessor.sequence) == 0:
                    bfb_distance = new_cost
                    start_backward = predecessor

            if cost[predecessor] == new_cost:
                # Equal sources at different indices yield parallel edges;
                # the first index is kept.
                G[predecessor].setdefault(node, operation)

    return G, start_backward, bfb_distance


def BFB_dist(sources, target):
    """Infer every history over ``sources`` reaching ``target`` that is
    minimum first in the number of t-BFB cycles and then in the total
    number of cycles.

    The returned :class:`~bfbdist.bfb_graph.BFBGraph` contains only edges
    between consecutive distance layers. Its distance is the number of
    cycles, counting the t-BFB cycles, zero when the target is empty, and
    ``math.inf`` when the target is unreachable. A source may be appended
    any number of times. The search runs backward from the target and
    deepens the bound on reversed t-BFB cycles from one upward; the first
    bound reaching the empty sequence is the minimum number of t-BFB
    cycles, and every history found under it uses exactly that many.
    """
    sources = validate_sources(sources)
    start = start_state()
    if len(target) == 0:
        return BFBGraph({start: {}}, sources, target, 0)

    for max_translocations in range(1, len(target) + 1):
        G, start_backward, bfb_distance = _bounded_search(
            sources, target, max_translocations
        )
        if start_backward is not None:
            adjacency = _relabel_forward(
                G, start_backward, start_backward.translocations
            )
            return BFBGraph(adjacency, sources, target, bfb_distance)
    return BFBGraph({start: {}}, sources, target, math.inf)
