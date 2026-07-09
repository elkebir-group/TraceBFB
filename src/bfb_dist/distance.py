from collections import defaultdict, deque
import math

def get_paths(G, source):
    if source not in G or not G[source]:
        return [[source]]
    return [[source] + p for t in G[source] for p in get_paths(G, t)]

def BFBDist(source, target):
    D = deque([target])
    G = defaultdict(set)
    cost = {target: 0}
    min_cost = math.inf
    
    while D:
        t = D.popleft()
        if cost[t] > min_cost:
            break
        for s in t.previous_bfb():
            new_cost = cost[t] + 1
            if s not in cost:
                cost[s] = new_cost
                D.append(s)
                G[s].add(t)
                if s == source:
                    min_cost = new_cost
            elif cost[s] == new_cost:
                G[s].add(t)
                
    return G, cost, min_cost, get_paths(G, source)