from collections import defaultdict, deque
import math
import random
import os
import networkx as nx


CHROMOSOME_LENGTHS = {
    "chr1": 248387328,
    "chr2": 242696752,
    "chr3": 201105948,
    "chr4": 193574945,
    "chr5": 182045439,
    "chr6": 172126628,
    "chr7": 160567428,
    "chr8": 146259331,
    "chr9": 150617247,
    "chr10": 134758134,
    "chr11": 135127769,
    "chr12": 133324548,
    "chr13": 113566686,
    "chr14": 101161492,
    "chr15": 99753195,
    "chr16": 96330374,
    "chr17": 84276897,
    "chr18": 80542538,
    "chr19": 61707364,
    "chr20": 66210255,
    "chr21": 45090682,
    "chr22": 51324926,
    "chrX": 154259566,
    "chrY": 62460029,
    "chrM": 16569,
}

def prRed(s): return f"\033[91m{s}\033[00m"
def prGreen(s): return f"\033[92m{s}\033[00m"
def prYellow(s): print("\033[93m {}\033[00m".format(s))
def prLightPurple(s): print("\033[94m {}\033[00m".format(s))
def prPurple(s): return f"\033[95m{s}\033[00m"
def prCyan(s): return f"\033[96m {s}\033[00m"
def prLightGray(s): print("\033[97m {}\033[00m".format(s))
def prBlack(s): print("\033[90m {}\033[00m".format(s))


class Segment:

    def __init__(self, chromosome, start, end, complement):
        self.chr = chromosome
        self.start = start
        self.end = end
        self.complement = complement

    def __eq__(self, other):
        return self.chr == other.chr and self.start == other.start and \
        self.end == other.end and self.complement == other.complement

    def __str__(self):
        return f'{self.chr}:{self.start}-{self.end}[R]' if not self.complement else f'{self.chr}:{self.start}-{self.end}[C]'

    __repr__ = __str__

    def reverse_complement(self):
        return Segment(self.chr, self.end, self.start, not self.complement)



class LetterAsSegment:

    def __init__(letter):
        self.letter = letter

    def reverse_complement(self):
        return letter.swapcase()



class Sequence:

    def __init__(self, segments=[]):
        self.segments = segments

    def __eq__(self, other):
        if len(self) == len(other):
            for x, y in zip(self.segments, other.segments):
                if x != y:
                    return False
        return True

    def __str__(self):
        return ';'.join([str(segment) for segment in self.segments])

    def __getitem__(self, i):
        if isinstance(i, slice):
            return Sequence(self.segments[i])
        else:
            self.segments[i]

    def __len__(self):
        return len(self.segments)

    def __add__(a, b):
        if isinstance(b, Sequence):
            return Sequence(a.segments + b.segments)
        elif isinstance(b, Segment):
            return Sequence(a.segments + [b])

    __repr__ = __str__

    def __hash__(self):
        return hash(str(self))

    def reverse_complement(self):
        return Sequence([segment.reverse_complement() for segment in self.segments[::-1]])

    def is_prefix_of(superseq):
        for i in range(len(self.segments)):
            if self.segments[i] != superseq[i]:
                return False
        return True

    def previous_bfb(self):
        pot_sources = set()
        l = len(self.segments)
        for i in range((l + 1) // 2, l):
            for x in range(0, 2 * i - l + 1):
                y = x + l - i
                if self[x:y] == self[i:].reverse_complement():
                    pot_sources.add(self[:i])
                    pot_sources.add(self[:y]+self[y:i].reverse_complement())
                    # print(self[:x], '------',self[x:y], '------',self[y:i], '------',self[i:], '->', self[i:].reverse_complement())
        return pot_sources

def simulate_reference_sequence(n_segments, min_segment_len=1000, max_segment_len=5000, starting_position=None):
    chromosome = random.choice(list(CHROMOSOME_LENGTHS.keys()))
    if starting_position is None:
        x = random.randint(0, CHROMOSOME_LENGTHS[chromosome] - n_segments * max_segment_len)
    else:
        x = starting_position
    ref = Sequence()
    for cycle in range(n_segments):
        y = x + random.randint(min_segment_len, max_segment_len)
        ref += Segment(chromosome, x, y, False)
        x = y
    return ref

def simulate_bfb(reference, num_cycles, verbose=True):
    template = reference
    ops = [reference]

    num_cycle = 0
    while num_cycle < num_cycles:
        x = random.randint(0, len(template) - 1)
        y = random.randint(x + 1,len(template) )
        num_cycle += 1
        p, q, r = template[:x], template[x:y], template[y:]
        if random.random() > 0.5:
            template = template + q.reverse_complement()
            ops.append(template)
        else:
            template = p + q + r.reverse_complement() + q.reverse_complement()
            ops.append(template)
        if verbose:
            print(prRed(p)+';'+prGreen(q)+';'+prPurple(r))
    if verbose:
        print(template)
    return ref, template, ops


def simulate_nested_bfb(reference, num_cycles, verbose=True):
    template = reference
    ops = [reference]

    num_cycle = 0
    while num_cycle < num_cycles:
        x = random.randint(0, len(template) - 1)
        y = random.randint(x + 1,len(template) )
        num_cycle += 1
        p, q, r = template[:x], template[x:y], template[y:]
        if random.random() > 0.5:
            template = template + q.reverse_complement()
            ops.append(template)
        else:
            template = p + q + r.reverse_complement() + q.reverse_complement()
            ops.append(template)
        if verbose:
            print(prRed(p)+';'+prGreen(q)+';'+prPurple(r))
    if verbose:
        print(template)
    return ref, template, ops

def get_paths(G, source):
    if source not in G:
        return [[source]]
    paths = []
    for t in G[source]:
        paths += [[source] + p for p in get_paths(G, t)]
    print(source, len(paths))
    return paths

def BFBDist(source, target):
    D = deque()
    G = defaultdict(set)
    cost = {target: 0}
    D.append(target)
    min_cost = math.inf
    while len(D) > 0:
        t = D.popleft()
        if cost[t] > min_cost:
            break
        for s in t.previous_bfb():
            if s not in cost or cost[s] == cost[t] + 1:
                if s not in cost:
                    D.append(s)
                G[s].add(t)
                cost[s] = cost[t] + 1
                if s == source:
                    min_cost = cost[s]
    return G, cost, min_cost, get_paths(G, source)