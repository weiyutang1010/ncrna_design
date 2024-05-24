import os, sys
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

FOLDER_PATH = "../results/sampling_pyx_targeted_sm_softmax_eps_075_adam_kmers"
DATA_PATH = "../data/eterna/eterna100.txt"
# puzzle_ids = [33, 47, 14]
# puzzle_ids = [77, 35]
# puzzle_ids = [69, 9]
# puzzle_ids = [43, 77]
# puzzle_ids = [23]
# puzzle_ids = [14, 33, 43, 47]
# puzzle_ids = [14, 33, 43, 47, 77, 69]

with open(DATA_PATH, 'r') as f:
    lines = f.read().split('\n')
    ids = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
    lengths = [int(len(line.split(' ')[1])) for line in lines if len(line) > 0]

def seqs_plot(seqs, avg_entropy, p_id, n):
    lines = []
    fig, ax1 = plt.subplots()
    for k in seqs:
        x = [i+1 for i in range(len(seqs[k]))]
        line, = ax1.plot(x, seqs[k], label=f'{k}-mers')
        lines.append(line)

    ax1.set_ylabel('Number of Unique Sequences')
    ax1.set_xlabel('Steps')
    
    # Average Positional Entropy
    x = [i+1 for i in range(len(avg_entropy))]
    ax2 = ax1.twinx()
    line, = ax2.plot(x, avg_entropy, linestyle='--', color='black', label=f'Pos Entropy')
    lines.append(line)
    ax2.set_ylabel('Average Positional Entropy')
    
    # sorted legends
    labels = [line.get_label() for line in lines]
    final_values = [line.get_ydata()[-1] for line in lines]
    sorted_lines_labels = sorted(zip(final_values, lines, labels), key=lambda x: -x[0])
    sorted_final_values, sorted_lines, sorted_labels = zip(*sorted_lines_labels)
    plt.legend(sorted_lines, sorted_labels, loc='upper left')

    plt.title(f'Puzzle {p_id} (n = {n})')

    save_path = f"./unique_seqs/{p_id}.png"
    plt.savefig(save_path, dpi=400, bbox_inches='tight')
    print(f"kmers unique seqs plot saved to {save_path}", file=sys.stderr)
    plt.close()

def ratio_plot(ratio, avg_entropy, p_id, n):
    lines = []
    fig, ax1 = plt.subplots()
    for k in ratio:
        x = [i+1 for i in range(len(ratio[k]))]
        line, = ax1.plot(x, ratio[k], label=f'{k}-mers')
        lines.append(line)

    # ax1.legend()
    ax1.set_ylabel('Uniqueness Ratio')
    ax1.set_xlabel('Steps')
    
    # Average Positional Entropy
    x = [i+1 for i in range(len(avg_entropy))]
    ax2 = ax1.twinx()
    line, = ax2.plot(x, avg_entropy, linestyle='--', color='black', label=f'Pos Entropy')
    lines.append(line)
    ax2.set_ylabel('Average Positional Entropy')
    
    # sorted legends
    labels = [line.get_label() for line in lines]
    final_values = [line.get_ydata()[-1] for line in lines]
    sorted_lines_labels = sorted(zip(final_values, lines, labels), key=lambda x: -x[0])
    sorted_final_values, sorted_lines, sorted_labels = zip(*sorted_lines_labels)
    plt.legend(sorted_lines, sorted_labels, loc='upper left')

    plt.title(f'Puzzle {p_id} (n = {n})')

    save_path = f"./unique_ratio/{p_id}.png"
    plt.savefig(save_path, dpi=400, bbox_inches='tight')
    print(f"kmers uniqueness ratio plot saved to {save_path}", file=sys.stderr)
    plt.close()

def uniqueness_plot(step_ratios):

    max_len = len(step_ratios[2000])
    print("Num Puzzles: ", max_len)
    print("Last Length: ", lengths[max_len - 1])
    for step, ratio in step_ratios.items():
        x = [i for i in range(max_len)]
        plt.bar(x, ratio[:max_len], label=f'Step {step}')
        print(f"Step: {step}, Average Uniqueness Ratio: {np.mean(ratio)}")
        # plt.plot(x, ratio, linestyle='', marker='x', label=f'Step {step}')

    plt.title('Uniqueness Ratio of Eterna100 Puzzles')
    plt.xlabel('Puzzle Lengths')
    plt.ylabel('Uniqueness Ratio (full sequence)')
    xtickslabels = [lengths[i] for i in np.arange(0, max_len, 5)]
    plt.xticks(np.arange(0, max_len, 5), xtickslabels)
    plt.legend()

    save_path = f"./uniqueness.png"
    plt.savefig(save_path, dpi=400, bbox_inches='tight')
    print(f"Puzzles uniqueness ratio plot saved to {save_path}", file=sys.stderr)
    plt.close()

def parse_file(puzzle_ids):
    step_ratios = {500: [], 1000: [], 1500: [], 2000: []}

    for p_id in puzzle_ids:
        file_path = f'{FOLDER_PATH}/{p_id}.txt'

        if not os.path.exists(file_path):
            continue

        k_arr = []
        ratio = {}
        seqs = {}

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

                    if k in seqs:
                        seqs[k].append(uniq_kmers)
                        ratio[k].append(uniq_kmers / total_count)
                    else:
                        seqs[k] = [uniq_kmers]
                        ratio[k] = [uniq_kmers / total_count]

                    if k == n and len(ratio[k]) in [500, 1000, 1500, 2000]:
                        step_ratios[len(ratio[k])].append(uniq_kmers / total_count)

                if line.startswith('entropy: '):
                    avg_entropy.append(float(line.split(': ')[1]))

        # seqs_plot(seqs, avg_entropy, p_id, n)
        # ratio_plot(ratio, avg_entropy, p_id, n)

    plt.figure(figsize=(12, 6))
    uniqueness_plot(step_ratios)

parse_file(ids)