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
    validate_sources,
)


def _relabel_forward(G, starts):
    """Return forward adjacency rows reachable from the empty-sequence states.

    The backward search labels a state with the set ``R`` of sources not yet
    stripped, which equals the sources appended so far plus the sources the
    history never uses. Each empty-sequence state in ``starts`` carries
    exactly the unused set ``N`` of the histories through it, so removing
    ``N`` from every label along those histories yields the forward
    convention in which the empty sequence has no source appended and every
    label is the set of sources appended so far.
    """
    adjacency = {}
    for start_node in starts:
        unused = start_node.used
        seen = {start_node}
        pending = [start_node]
        while pending:
            node = pending.pop()
            children = G.get(node)
            if not children:
                continue
            row = adjacency.setdefault(
                State(node.sequence, node.used & ~unused), {}
            )
            for child, operation in children.items():
                row[State(child.sequence, child.used & ~unused)] = operation
                if child not in seen:
                    seen.add(child)
                    pending.append(child)
    return adjacency or {start_state(): {}}


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
    target is empty, and ``math.inf`` when the target is unreachable. Each
    source is appended at most once; sources a shortest history does not
    need are left unused.
    """
    sources = validate_sources(sources)
    start = start_state()
    if len(target) == 0:
        return BFBGraph({start: {}}, sources, target, 0)

    # The search runs backward from the target with every source available
    # and stops at the first layer holding an empty sequence; the sources
    # still available there are the ones the history leaves unused.
    full = end_state(sources, target)
    pending = deque([full])
    G = defaultdict(dict)
    cost = {full: 0}
    bfb_distance = math.inf
    starts = []

    while pending:
        node = pending.popleft()
        if cost[node] >= bfb_distance:
            break

        new_cost = cost[node] + 1
        for predecessor, operation in _predecessors(node, sources):
            if predecessor not in cost:
                cost[predecessor] = new_cost
                pending.append(predecessor)
                if len(predecessor.sequence) == 0:
                    bfb_distance = new_cost
                    starts.append(predecessor)

            if cost[predecessor] == new_cost:
                G[predecessor][node] = operation

    adjacency = _relabel_forward(G, starts) if starts else {start: {}}
    return BFBGraph(adjacency, sources, target, bfb_distance)
