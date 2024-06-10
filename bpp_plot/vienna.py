"""Given a sequence and structure, output its p(y | x), NED, is_mfe, is_umfe"""
"""Reference: https://github.com/shanry/SAMFEO/blob/main/utils/vienna.py"""

import os, sys
import numpy as np
import RNA

def extract_pairs(ss):
    pairs = list(range(len(ss)))
    stack = []
    for i, c in enumerate(ss):
        if c=='.':
            pass
        elif c=="(":
            stack.append(i)
        elif c==")":
            j = stack.pop()
            pairs[j] = i
            pairs[i] = j
        else:
            raise ValueError(f"wrong structure at position {i}: {c}")
    return pairs

def base_pair_probs(seq, sym=False, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = np.array(fc.bpp())[1:, 1:]
    if sym:
        bpp += bpp.T
        unpair = 1 - np.sum(bpp, axis=1)
        bpp[range(len(bpp)), range(len(bpp))] = unpair
    return bpp

def mfe(seq):
    fc = RNA.fold_compound(seq)
    ss = fc.mfe()
    return ss

# Print a subopt result as FASTA record
def print_subopt_result(structure, energy, data):
    ss_list = []
    if not structure == None:
        data['ss_list'].append((energy, structure))
        data['counter'] = data['counter'] + 1

def subopt(seq, e=0):
    subopt_data = { 'counter' : 0, 'sequence' : seq, 'ss_list': []}
    fc = RNA.fold_compound(seq)
    fc.subopt_cb(e, print_subopt_result, subopt_data)
    subopt_data['ss_list'] = sorted(subopt_data['ss_list'])
    return subopt_data

def prob(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    pr = fc.pr_structure(ss)
    return pr

def ensemble_defect(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    fc.bpp()
    ed = fc.ensemble_defect(ss)
    return ed

def position_defect(seq, ss, scale=True):
    bpp = base_pair_probs(seq, sym=True, scale=True)
    pairs = extract_pairs(ss)
    defect_pos = [1 - bpp[i, j] for i, j in enumerate(pairs)]

    return bpp, defect_pos
