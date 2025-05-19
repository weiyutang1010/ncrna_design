#!/usr/bin/env python3

import pandas as pd
import csv
import numpy as np
import RNA

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

def mfe(seq):
    fc = RNA.fold_compound(seq)
    ss = fc.mfe()
    return ss

def structural_dist(seq, ss):
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

def e_diff(seq, ss, ss_mfe=""):
    if ss_mfe == "":
        ss_mfe = mfe(seq)[0]
    return abs(energy(seq, ss_mfe) - energy(seq, ss))

undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]
unknown_ids = [68, 97, 100]

file = 'ddg_softmax'

# Read the CSV file into a DataFrame
df = pd.read_csv(f'old_data/{file}.csv')
df.rename(columns={'p(y | x)': 'p(y|x)'}, inplace=True)

print(df.columns)

designability = ['designable'] * 100
ids = df['ID']

for idx, pid in enumerate(ids):
    if pid in undesignable_ids:
        designability[idx] = 'undesignable'
    if pid in unknown_ids:
        designability[idx] = 'unknown'

df['Designability'] = designability

df = df.reindex(columns=['ID', 'Name', 'Length', 'Structure', 'Unpaired', 'Pairs', 'Total', 'Designability',
       'p(y|x)', 'p(y|x) seq', 'NED', 'NED seq', 'dist', 'dist seq', 'DDG',
       'DDG seq', 'is_mfe', 'mfe_seq', 'is_umfe', 'umfe_seq'])

del df['Total']

pyx_arr = []
ned_arr = []
dist_arr = []
ddg_arr = []

structs = df['Structure']
seqs = df['p(y|x) seq']
for seq, struct in zip(seqs, structs):
    pyx = prob(seq, struct)

    if pyx < 0.001:
        pyx_arr.append(f'{pyx:.1e}')
    else:
        pyx_arr.append(f'{pyx:.3f}') 

seqs = df['NED seq']
for seq, struct in zip(seqs, structs):
    ned = ensemble_defect(seq, struct)
    ned_arr.append(f'{ned:.3f}')

seqs = df['DDG seq']
for seq, struct in zip(seqs, structs):
    ddg = e_diff(seq, struct)
    ddg_arr.append(f'{ddg:.2f}')

df['p(y|x)'] = pyx_arr
df['NED'] = ned_arr
df['DDG'] = ddg_arr

df.to_csv(f'{file}.csv', index=False)

print("Done")
# print(df['p(y|x)'])