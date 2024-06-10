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

    return defect_pos

def mfe(seq):
    fc = RNA.fold_compound(seq)
    ss = fc.mfe()
    return ss

def structural_dist(seq, ss, ss_mfe):
    ss_mfe = mfe(seq)[0]
    stk = []
    mp = {}

    for j, c in enumerate(ss):
        if c == '(':
            stk.append(j)
        elif c == ')':
            i = stk.pop()
            mp[j] = i
            mp[i] = j
        else:
            mp[j] = -1

    dist = len(ss)
    for j, c in enumerate(ss_mfe):
        if c == '(':
            stk.append(j)
        elif c == ')':
            i = stk.pop()
            
            if mp[j] == i:
                dist -= 2
        else:
            if mp[j] == -1:
                dist -= 1

    return dist

def energy(seq, ss):
    fc = RNA.fold_compound(seq)
    return fc.eval_structure(ss)

def e_diff(seq, ss, ss_mfe):
    return abs(energy(seq, ss_mfe) - energy(seq, ss))

if __name__ == '__main__':
    # with open('input.txt', 'r') as f:
    #     lines = f.read().split('\n')
    #     n = len(lines)
        
    # seqs    = lines[:n//2]
    # structs = lines[n//2:]

    seq =    "GGACGGAACCGCGACGAAAAUGGACCGCGAAAGCAAAAGGGACACGAAGCAAAAGGGACCGCGACCAACGGACGGCAACCAAAAUCGACCGCGAGGAAAACCAAAAGGCCAAGCAAAAGGGAGAGGAGGCGCGAGCAAGGGAAACCGGGGGAGCAGGGAAAAGGGACCGCGAAAGCAAAAGGGACCCCAAGCAAAAUCAACCGCGAGCAACGGACGGGGACUAAAAGGGACAGCAAGGAAAACCGAAAGGCCGAGCAAAAUGGACCAGGACCGGGAGCAAGGAC"
    # seq = "GGGCAAAGCCCGGCACUGGAAACAGAGGCAAGGGACCCGAAAGGGGAGCGCCCAAAGGGCAACGGGAAACCCGGCUCCGGGAAACCGUCCCAAGAGAGGCGAAAGCGCGCCGCCCUAAGGGCAAGGUCAAAGACCGGCGCUGGAAACACUCUCAACCUCCUGGAAACACUGCCGCCCAAAGGGCUACGGCAAAGCCGGGCAGGCGAAAGCGGAGGAACCGGCGCGAAAGCGGGCCGCGGAAUCCGCAACUCCAAAGGAGGGCCCAGGAAACUGCCGGAAGGAGCCCGAAAGGGAGGCGCCCCAAGGGCAAGGGCAAAGCCCGCCUCGGGAAACCGCUCCAAGCCUCGCGAAAGCGUGCCCCGGAAACCGG"
    struct = "(..(((..((((..((....((..((((....))....))..))))..))....))..))((..((..((..((((..((....((..((((..((....((....))))..))....))..))))..))((..((..((....))((..((((..((....((..((((....))....))..))))..))....))..))))..))..))..))((..((....((..((((..((....((....))))..))....))..))))..))))..))..)..)"

    seqs = [seq]
    structs = [struct]

    for seq, struct in zip(seqs, structs):
        if seq == "" or struct == "":
            exit(0)

        if len(seq) != len(struct):
            print("Length doesn't match!")

        ss_mfe = mfe(seq)[0]
        subopt_data = subopt(seq)

        mfe_structs = [st for e, st in subopt_data['ss_list']]
        is_mfe = struct in mfe_structs
        is_umfe = is_mfe and subopt_data['counter'] == 1
        pyx = prob(seq, struct)
        ned = ensemble_defect(seq, struct)
        pos_def = position_defect(seq, struct)
        dist = structural_dist(seq, struct, ss_mfe)
        delta_delta_G = e_diff(seq, struct, ss_mfe)

        print(f"seq        : {seq}")
        print(f"structure  : {struct}")
        # print(f"mfe(seq) : {mfe_structs}")
        print(f"is_mfe     : {is_mfe}")
        print(f"is_umfe    : {is_umfe}")
        print(f"dist       : {dist}")
        print(f"ediff      : {delta_delta_G}")
        print(f"p(y | x)   : {pyx}")
        print(f"ned        : {ned}")
        # print(f"pos defect : {pos_def}")
        print(f"mfe_struct : {ss_mfe}")
        print()