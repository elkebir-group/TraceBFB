"""Nested and strongly-nested deletion-aware BFB distance inference."""

from collections import defaultdict, deque
import math

from bfb_dist.bfb_graph import BFBGraph


def _is_foldable(node, target):
    """Return whether node's tail beyond its LCP with target is foldable."""
    tail = node[node.lcp_length(target):]
    return len(tail) == 0 or tail.reverse_complement().is_infix_of(target)


def _truncate_G(G, source, target):
    """Return adjacency rows for nodes with a recorded path to ``target``.

    ``_search`` explores forward from ``source``, so every key of ``G`` is
    already reachable from ``source`` by construction; the missing property
    is reachability *to* ``target``. This walks backward from ``target``
    over ``G``'s reversed edges to find it -- the mirror image of
    ``distance._truncate_G``, which walks forward from ``source`` because
    ``BFB_dist`` explores backward from ``target`` instead. Not
    interchangeable with that helper despite the identical name/shape.
    """
    parents_of = defaultdict(list)
    for node, children in G.items():
        for child in children:
            parents_of[child].append(node)

    reachable = {target}
    pending = deque([target])
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
    return adjacency or {source: {}}


def _search(source, target, star):
    """Infer every minimum-length (strongly) nested BFB path source-to-target."""
    if source == target:
        return BFBGraph({source: {}}, source, target, 0)

    pending = deque([source])
    G = defaultdict(dict)
    cost = {source: 0}
    bfb_distance = math.inf

    while pending:
        node = pending.popleft()
        if cost[node] >= bfb_distance:
            break

        node_phi = node.lcp_length(target)
        new_cost = cost[node] + 1
        for child, x, y, is_forward_deletion in node.successors():
            if child.lcp_length(target) <= node_phi:
                continue
            if star and not _is_foldable(child, target):
                continue

            if child not in cost:
                cost[child] = new_cost
                pending.append(child)
                if child == target:
                    bfb_distance = new_cost

            if cost[child] == new_cost:
                G[node][child] = (x, y, is_forward_deletion)

    adjacency = _truncate_G(G, source, target)
    return BFBGraph(adjacency, source, target, bfb_distance)


def NBFB_dist(source, target):
    """Infer every minimum-length nested BFB path from source to target."""
    return _search(source, target, star=False)


def NBFB_star_dist(source, target):
    """Infer every minimum-length strongly nested BFB path from source to target."""
    return _search(source, target, star=True)
