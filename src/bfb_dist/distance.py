"""Minimum-distance BFB graph inference."""

from collections import defaultdict, deque
import math

from bfb_dist.bfb_graph import BFBGraph


def _truncate_G(G, source):
    """Return reachable adjacency rows, omitting terminal nodes."""
    adjacency = {}
    visited = {source}
    pending = deque([source])

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

    return adjacency or {source: {}}


def BFB_dist(source, target):
    """Infer every minimum-distance BFB path from source to target.

    The returned :class:`~bfb_dist.bfb_graph.BFBGraph` contains only edges
    between consecutive reverse-distance layers. Its distance is zero when
    source equals target and ``math.inf`` when target is unreachable.
    """
    if source == target:
        return BFBGraph({source: {}}, source, target, 0)

    pending = deque([target])
    G = defaultdict(dict)
    cost = {target: 0}
    bfb_distance = math.inf

    while pending:
        node = pending.popleft()
        if cost[node] >= bfb_distance:
            break

        new_cost = cost[node] + 1
        for predecessor, x, y, is_forward_deletion in node.previous_bfb():
            if y == len(predecessor):
                is_forward_deletion = False

            if predecessor not in cost:
                cost[predecessor] = new_cost
                pending.append(predecessor)
                if predecessor == source:
                    bfb_distance = new_cost

            if cost[predecessor] == new_cost:
                G[predecessor][node] = (
                    x,
                    y,
                    is_forward_deletion,
                )

    adjacency = _truncate_G(G, source)
    return BFBGraph(adjacency, source, target, bfb_distance)
