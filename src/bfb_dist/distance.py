"""Deletion-aware BFB distance over a set of sources (backward search)."""

from collections import defaultdict, deque
import math

from bfb_dist.bfb_graph import BFBGraph
from bfb_dist.history import (
    Operation,
    State,
    end_state,
    mask_indices,
    start_state,
    total_source_length,
    validate_sources,
)


def _truncate_G(G, start):
    """Return reachable adjacency rows, omitting terminal nodes."""
    adjacency = {}
    visited = {start}
    pending = deque([start])

    while pending:
        node = pending.popleft()
        children = G.get(node)
        if not children:
            continue

        adjacency[node] = children
        for child in children:
            if child not in visited:
                visited.add(child)
                pending.append(child)

    return adjacency or {start: {}}


def _predecessors(node, sources):
    """Yield ``(predecessor, operation)`` pairs producing ``node``.

    Reversed p-, i-, and f-BFB cycles keep the set of appended sources.
    A reversed t-BFB cycle strips a source that is a suffix of the sequence
    and removes the source from the set.
    """
    sequence = node.sequence
    for predecessor, x, y, is_forward_deletion in sequence.previous_bfb():
        operation = Operation.fold_back(
            x, y, is_forward_deletion, len(predecessor)
        )
        yield State(predecessor, node.used), operation

    for index in mask_indices(node.used):
        predecessor = sequence.strip_suffix(sources[index])
        if predecessor is not None:
            yield (
                State(predecessor, node.used & ~(1 << index)),
                Operation.translocation(index),
            )


def BFB_dist(sources, target):
    """Infer every minimum-length history over ``sources`` reaching ``target``.

    The returned :class:`~bfb_dist.bfb_graph.BFBGraph` contains only edges
    between consecutive distance layers. Its distance is the number of
    cycles of a shortest history, counting the t-BFB cycles, zero when the
    target is empty and there is no source, and ``math.inf`` when the target
    is unreachable.
    """
    sources = validate_sources(sources)
    start, end = start_state(), end_state(sources, target)
    if start == end:
        return BFBGraph({start: {}}, sources, target, 0)
    if total_source_length(sources) > len(target):
        return BFBGraph({start: {}}, sources, target, math.inf)

    pending = deque([end])
    G = defaultdict(dict)
    cost = {end: 0}
    bfb_distance = math.inf

    while pending:
        node = pending.popleft()
        if cost[node] >= bfb_distance:
            break

        new_cost = cost[node] + 1
        for predecessor, operation in _predecessors(node, sources):
            if predecessor not in cost:
                cost[predecessor] = new_cost
                pending.append(predecessor)
                if predecessor == start:
                    bfb_distance = new_cost

            if cost[predecessor] == new_cost:
                G[predecessor][node] = operation

    adjacency = _truncate_G(G, start)
    return BFBGraph(adjacency, sources, target, bfb_distance)
