import os, sys
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

FOLDER_PATH = "../results/sampling_pyx_targeted_sm_softmax_eps_075_adam_kmers"

puzzle_ids = [33, 47, 14]
# puzzle_ids = [77, 35]
# puzzle_ids = [69, 9]
# puzzle_ids = [43, 77]
# puzzle_ids = [14, 33, 43, 47]
# puzzle_ids = [14, 33, 43, 47, 77, 35]

for p_id in puzzle_ids:
    file_path = f'{FOLDER_PATH}/{p_id}.txt'

    k_arr = []
    ratio = {}
    n = 0

    avg_entropy = []
    with open(file_path, 'r') as f:
        lines = f.read().split('\n')
        n = len(lines[0])

        for line in lines:
            if line.startswith('k: '):
                line_split = line.split(', ')

                k = int(line_split[0].split(': ')[1])
                uniq_kmers = int(line_split[1].split(': ')[1])
                total_count = int(line_split[2].split(': ')[1])

                if k in ratio:
                    ratio[k].append(uniq_kmers / total_count)
                else:
                    ratio[k] = [uniq_kmers / total_count]

            if line.startswith('entropy: '):
                avg_entropy.append(float(line.split(': ')[1]))

    lines = []
    fig, ax1 = plt.subplots()
    for k in reversed(ratio):
        if k != 3:
            x = [i+1 for i in range(len(ratio[k]))]
            line, = ax1.plot(x, ratio[k], label=f'{k}-mers')
            lines.append(line)

    # ax1.legend()
    ax1.set_ylabel('Unique Ratio')
    ax1.set_xlabel('Steps')
    
    x = [i+1 for i in range(len(avg_entropy))]
    ax2 = ax1.twinx()
    line, = ax2.plot(x, avg_entropy, linestyle='--', color='black', label=f'Average Positional Entropy')
    lines.append(line)
    ax2.set_ylabel('Average Positional Entropy')
    
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels)

    plt.title(f'Puzzle {p_id} (n = {n})')

    save_path = f"./plots/{p_id}.png"
    plt.savefig(save_path, dpi=400, bbox_inches='tight')
    print(f"kmers plot saved to {save_path}", file=sys.stderr)
    plt.close()