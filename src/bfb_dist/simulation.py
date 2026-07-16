import random
from bfb_dist.utils import CHROMOSOME_LENGTHS, format_red, format_green, format_purple
from bfb_dist.sequence import Segment, Sequence

def _set_seed(seed):
    random.seed(seed)

def simulate_reference_sequence(n_segments, min_segment_len=1000, max_segment_len=5000, starting_position=None):
    chromosome = random.choice(list(CHROMOSOME_LENGTHS.keys()))
    if starting_position is None:
        max_bound = CHROMOSOME_LENGTHS[chromosome] - n_segments * max_segment_len
        x = random.randint(0, max(0, max_bound))
    else:
        x = starting_position
    ref = Sequence()
    for _ in range(n_segments):
        y = x + random.randint(min_segment_len, max_segment_len)
        ref += Segment(chromosome, x, y, False)
        x = y
    return ref

def simulate_bfb(reference, num_cycles, verbose=True):
    template = reference
    ops = [reference]
    for n_cycle in range(num_cycles):
        l = len(template)
        x = random.randint(0, l - 1) if l > 1 else 0
        y = random.randint(x + 1, l) if l > 0 else l
        p, q, r = template[:x], template[x:y], template[y:]
        
        if random.random() > 0.5:
            template = template + q.reverse_complement()
        else:
            template = p + q + r.reverse_complement() + q.reverse_complement()
        ops.append(template)
        
        if verbose:
            print('R: ' if n_cycle == 0 else '-> ', ';'.join(s for s in [format_red(p), format_green(q), format_purple(r)] if len(s) > 16))
    if verbose:
        print('-> ', template)
    return reference, template, ops