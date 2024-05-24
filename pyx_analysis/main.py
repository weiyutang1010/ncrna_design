import os, sys
import numpy as np
import matplotlib.pyplot as plt

import RNA

puzzles = []

with open('../data/eterna/eterna100.txt', 'r') as data:
    lines = data.read().split('\n')
    puzzles = [line.split(' ')[1] for line in lines if len(line) > 0]

sampling = []
samfeo = []
samfeo_nosm = []

def get_boltz_prob(seq, ss, scale=False):
    """viennaRNA boltzmann probability"""
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    Q = fc.pf()
    pr = fc.pr_structure(ss)
    return pr

def parse(file_path):
    res = []
    with open(file_path, 'r') as f:
        lines = f.read().split('\n')
        for idx, line in enumerate(lines):
            res.append(get_boltz_prob(line, puzzles[idx]))

    return res

sampling = parse("./sampling.txt")
samfeo = parse("./samfeo.txt")
samfeo_nosm = parse("./samfeo_nosm.txt")

# print(np.mean(sampling))
# print(np.mean(samfeo))
# print(np.mean(samfeo_nosm))

x = np.arange(10, 101, 10)

sampling_bar = []
samfeo_bar = []
samfeo_nosm_bar = []

prev = 0
for i in range(10, 101, 10):
    sampling_bar.append(np.mean(sampling[prev:i]))
    samfeo_bar.append(np.mean(samfeo[prev:i]))
    samfeo_nosm_bar.append(np.mean(samfeo_nosm[prev:i]))
    prev = i

width = 2.5

fig, ax = plt.subplots(figsize=(14, 6))
bars1 = ax.bar(x - width, sampling_bar, width, label='Sampling (Softmax, Adam)', color='lightblue', edgecolor='black')
bars2 = ax.bar(x, samfeo_bar, width, label='SAMFEO', color='lightgreen', edgecolor='black')
bars3 = ax.bar(x + width, samfeo_nosm_bar, width, label='SAMFEO (no sm)', color='lightcoral', edgecolor='black')

# Add labels and title
ax.set_xlabel('Puzzles')
ax.set_ylabel('Average p(y | x)')
ax.set_title('Average p(y | x) for every 10 puzzles')
ax.set_xticks(x)
plt.legend()

plt.savefig('./plot.png', dpi=400, bbox_inches='tight')
plt.close()