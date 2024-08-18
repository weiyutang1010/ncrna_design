#!/usr/bin/env python3

import os, sys
import numpy as np
import pandas as pd
from scipy.stats import gmean
from vienna import *
from collections import defaultdict
import concurrent.futures


import matplotlib.pyplot as plt

import RNA

files = [
    'eterna100_nemo_rp_0.csv',
    'eterna100_nemo_rp_1.csv',
    'eterna100_nemo_rp_2.csv',
    'eterna100_nemo_rp_3.csv',
    'eterna100_nemo_rp_4.csv',
]

class Result:
    def __init__(self):
        self.best_pyx, self.best_pyx_seq, self.pyx_step_found, self.pyx_file = -0.1, "", -1, ""
        self.best_ned, self.best_ned_seq = 1.1, ""
        self.best_dist, self.best_dist_seq = 1000, ""
        self.best_ddg, self.best_ddg_seq = 100000, ""
        self.mfe_found, self.mfe_seq = False, ""
        self.umfe_found, self.umfe_seq = False, ""

results = defaultdict(Result) # {id: (metrics*)}

def update_result(pid, seq, struct):
    seq, pyx, ned, is_mfe, is_umfe, dist, ddg = eval_seq(seq, struct)

    if pyx > results[pid].best_pyx:
        results[pid].best_pyx = pyx
        results[pid].best_pyx_seq = seq

    if ned < results[pid].best_ned:
        results[pid].best_ned = ned
        results[pid].best_ned_seq = seq

    if dist < results[pid].best_dist:
        results[pid].best_dist = dist
        results[pid].best_dist_seq = seq

    if ddg < results[pid].best_ddg:
        results[pid].best_ddg = ddg
        results[pid].best_ddg_seq = seq

    if is_mfe == True:
        results[pid].mfe_found = True
        results[pid].mfe_seq = seq

    if is_umfe == True:
        results[pid].umfe_found = True
        results[pid].umfe_seq = seq

results = defaultdict(Result)

for file in files:
    df = pd.read_csv(file)

    for i in range(100):
        seq = df['rna'].iloc[i]
        struct = df['structure'].iloc[i]

        update_result(i, seq, struct)


with open('../data/eterna/eterna100.txt', 'r') as data:
    lines = data.read().split('\n')
    ids = [int(line.split(' ')[0]) for line in lines if len(line) > 0]

probs = [result.best_pyx for result in results.values()]

temp = probs[:]
probs[95] = temp[98]
probs[96] = temp[95]
probs[97] = temp[99]
probs[98] = temp[97]
probs[99] = temp[96]

undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]
probs_no_undesignable = [prob for pid, prob in zip(ids, probs) if pid not in undesignable_ids]
print(gmean(probs_no_undesignable))

# statements = []
# for i in range(100):
    # is_mfe = "is_mfe" if results[i].mfe_found else ""
    # mfe_seq = results[i].mfe_seq

    # is_umfe = "is_umfe" if results[i].umfe_found else ""
    # umfe_seq = results[i].umfe_seq

    # line = f"{results[i].best_pyx},{results[i].best_pyx_seq},{results[i].best_ned},{results[i].best_ned_seq},{results[i].best_dist},{results[i].best_dist_seq},{results[i].best_ddg},{results[i].best_ddg_seq},{is_mfe},{mfe_seq},{is_umfe},{umfe_seq}"

    # statements.append(line)

# temp = statements[:]
# statements[95] = temp[98]
# statements[96] = temp[95]
# statements[97] = temp[99]
# statements[98] = temp[97]
# statements[99] = temp[96]

# for line in statements:
#     print(line)
