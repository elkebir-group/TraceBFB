from collections import defaultdict, deque
import math

def get_paths(G, source):
    if source not in G or not G[source]:
        return [[source]]
    return [[source] + p for t in G[source] for p in get_paths(G, t)]

def _truncate_G(G, source):
    nodes = set([source])
    D = deque([source])
    while D:
        s = D.popleft()
        for t in G[s]:
            if t not in nodes:
                nodes.add(t)
                D.append(t)
    return {i: set([j:k for (j,k) in G[i].items() if j in nodes]) for i in nodes}

def BFB_dist(source, target):
    D = deque([target])
    G = defaultdict(dict)
    cost = {target: 0}
    min_cost = math.inf
    
    while D:
        t = D.popleft()
        if cost[t] > min_cost:
            break
        for s, s_decomp in t.previous_bfb().items():
            new_cost = cost[t] + 1
            if s not in cost:
                cost[s] = new_cost
                D.append(s)
                G[s][t] = s_decomp
                if s == source:
                    min_cost = new_cost
            elif cost[s] == new_cost:
                G[s][t] = s_decomp
                
    return _truncate_G(G, source), min_cost