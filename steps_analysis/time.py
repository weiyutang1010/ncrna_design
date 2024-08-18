#!/usr/bin/env python3

import os, sys
import numpy as np
import matplotlib.pyplot as plt

with open('../data/eterna/eterna100.txt', 'r') as data:
    lines = data.read().split('\n')
    ids = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
    structs = [line.split(' ')[1] for line in lines if len(line) > 0]


time_sec = []

for idx, pid in enumerate(ids):
    with open(f'../results/sampling_time/{pid}.txt', 'r') as f:
        lines = f.read().split('\n')
        
        if len(lines) > 2 and lines[-2].startswith('Total Time: '):
            time = float(lines[-2].split(': ')[1])
            time_sec.append(time / 3600)
        else:
            break

x = [len(structs[i]) for i in range(len(time_sec))]
print(len(x))

plt.rcParams["figure.figsize"] = [14, 5.50]

plt.plot(x, np.array(time_sec), linestyle='', marker='o', markerfacecolor='None', ms=10, alpha=0.8)

# determine best fit line
z = np.polyfit(x, time_sec, 3)

p = np.poly1d(z)
xl = np.linspace(min(x), max(x), 100)
yl = [p(xx)  for xx in xl]

text = f'Line of Best Fit (Cubic)'
plt.plot(xl, yl, label=text)
plt.legend(fontsize=16)

# xticks = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# plt.xticks(xticks)
plt.xlim(0, 410)
plt.xlabel(r'Puzzle Length', fontsize=15)
plt.ylabel('Time (hour)', fontsize=15)
plt.tick_params(labelsize=14)
plt.savefig("time.pdf", bbox_inches='tight')