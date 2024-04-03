import os, sys
import numpy as np
import matplotlib.pyplot as plt

with open('../data/eterna/eterna_n50.txt', 'r') as data:
    lines = data.read().split('\n')
    puzzles = [(line.split(' ')[0], line.split(' ')[1]) for line in lines if len(line) > 0]

init = sys.argv[1]
folder = f'../results/boxplot_{init}'

marker_props = dict(marker='.', markerfacecolor='black', markersize=2, linestyle='none')

for p_id, struct in puzzles:
    mean_val = []
    results = []
    x_values = []
    with open(f'{folder}/{p_id}.txt', 'r') as f:
        lines = f.read().split('\n')

        start_idx = lines.index('best samples')

        step = 0
        for idx in range(start_idx, len(lines), 2503):
            if lines[idx] == 'best samples':
                if step % 200 == 0 or step == 1999:
                    j = idx + 1
                    temp = []
                    while len(lines[j]) > 0:
                        temp.append(float(lines[j].split(' ')[1]))
                        j += 1
                    results.append(temp[:])
                    mean_val.append(np.mean(temp))
                    x_values.append(step)
                step += 1

        plt.figure(figsize=(12, 5))  # Adjust width and height as needed
        plt.xlabel('Step')
        plt.ylabel('Probability')
        plt.title(f'Samples Probability, init={init}, id={p_id}')
        plt.boxplot(results, widths=60, positions=x_values, flierprops=marker_props)
        plt.scatter(x_values, mean_val, marker='.', s=50, color='red')
        plt.savefig(f'./boxplot_{init}_p{p_id}.png', format='png', bbox_inches='tight')
