import sys
import RNA
from itertools import product

import threading
import concurrent.futures

lock = threading.Lock()

def get_boltz_prob(seq, ss, scale=True):
    """viennaRNA boltzmann probability"""
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    pr = fc.pr_structure(ss)

    with lock:
        print(f"{seq} {pr}")

def valid_seq(y, targeted=False):
    x = []

    pairs = []
    unpairs = []
    stack = []
    for j, c in enumerate(y):
        if c == '(':
            stack.append(j)
        elif c == ')':
            pairs.append((stack.pop(), j))
        else:
            unpairs.append(j)

    if targeted:
        seq = ['A'] * len(y)
        for combination in product(['CG', 'GC'], repeat=len(pairs)):
            for idx, pair in enumerate(combination):
                i, j = pairs[idx]
                seq[i] = pair[0]
                seq[j] = pair[1]
            x.append("".join(seq))
    else:
        seq = ['A'] * len(y)
        for comb1 in product(['CG', 'GC', 'AU', 'UA', 'GU', 'UG'], repeat=len(pairs)):
            for comb2 in product(['A', 'C', 'G', 'U'], repeat=len(unpairs)):
                for idx, pair in enumerate(comb1):
                    i, j = pairs[idx]
                    seq[i] = pair[0]
                    seq[j] = pair[1]
                
                for idx, c in enumerate(comb2):
                    seq[unpairs[idx]] = c

                x.append("".join(seq))

    return x


rna_struct = sys.argv[1]
print(rna_struct)
with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
    futures = [executor.submit(get_boltz_prob, seq, rna_struct) for seq in valid_seq(rna_struct)]
    concurrent.futures.wait(futures)