import os, sys
import numpy as np
import matplotlib.pyplot as plt

with open("../data/eterna/eterna_n50.txt", 'r') as f:
    data = f.read().split('\n')
    data = [(int(line.split(' ')[0]), line.split(' ')[1]) for line in data if len(line) > 0]

steps = [x+1 for x in range(1500)]
for p_id, struct in data:
    with open(f"./p{p_id}.txt", 'r') as f:
        lines = f.read().split('\n')
        if len(lines) < 1500:
            continue
        n = len(struct)

        var = np.array([float(line.split(': ')[1]) for line in lines if len(line) > 0])
        plt.plot(steps, var, alpha=0.85, label=f'p{p_id}')

plt.title("Sample Variance vs. Steps (Uniform Initialization)")
plt.xlabel("Steps")
plt.ylabel("Sample Variance")
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(f'./variance.png', format='png', bbox_inches='tight')