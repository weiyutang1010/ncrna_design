import os, sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
from collections import defaultdict

import RNA

def get_boltz_prob(seq, ss, scale=True):
    """viennaRNA boltzmann probability"""
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

def get_best_sample(lines, rna_struct):
    samples = []
    for line in lines:
        seq = line.split(' ')[0]
        samples.append((seq, get_boltz_prob(seq, rna_struct)))

    samples.sort(key= lambda x: x[1])
    return samples[-1]

def process_result_file(rna_id, result_file, args):
    lines = result_file.read().split('\n')

    if len(lines) < 5:
        exit(0)

    n, rna_struct = len(lines[0]), lines[0]

    init = lines[1].split(', ')[1].split(': ')[1]

    learning_rate = lines[2].split(', ')[0].split(': ')[1]
    sample_size = lines[3].split(', ')[0].split(': ')[1]

    best_k = 3
    # best_k = int(lines[3].split(', ')[2].split(': ')[1])

    time = 0.0
    if lines[-2].startswith("Total Time: "):
        time = float(lines[-2].split(': ')[1])

    objs, seqs, pyx, best_samples = [], [], [], []
    approx_prob = []
    boxplot = []
    
    # File reading
    for idx, line in enumerate(lines):
        if line.startswith("Boxplot: "):
            values = line.split(': ')[1].split(' ')
            values = [float(value) for value in values if len(value) > 0]
            boxplot.append(values)

        if line.startswith("step:"):
            values = line.split(', ')
            objs.append(float(values[1].split(': ')[1]))
            approx_prob.append(float(values[2].split(': ')[1]))
            seqs.append(values[3].split(': ')[1])

        if line.startswith("best samples"):
            sample_lines, j = [], idx + 1
            for j in range(idx + 1, idx + 1 + best_k):
                sample_lines.append(lines[j])

            best_samples.append(get_best_sample(sample_lines[:2], rna_struct))


    if len(objs) < 1:
        exit(0)
    
    sample_probs = [sample[1] for sample in best_samples]
    objs_exp = [np.exp(-1 * obj) for obj in objs]

    best_seq, best_score = '', 0.
    prev_seq, prev_score = '', 0.
    for idx, seq in enumerate(seqs):
        if seq != prev_seq:
            prev_seq = seq
            prev_score = get_boltz_prob(seq, rna_struct)
            pyx.append(prev_score)
            if prev_score > best_score:
                best_seq = prev_seq
                best_score = prev_score
        else:
            pyx.append(prev_score)

    print("Best Integral Solution (seq, score): ", best_seq, best_score)
    best_sample_seq, best_sample_prob = max(best_samples, key=lambda x: x[1])
    print("Best Sampled Solution (seq, score): ", best_sample_seq, best_sample_prob)

    plt.rcParams["figure.figsize"] = [7.50, 4.50]
    plt.rcParams["figure.autolayout"] = True

    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Step')
    ax1.set_ylabel('Boltzmann Probability')

    # box plot
    num_steps = len(objs)
    x_values = [x for x in range(0, num_steps, num_steps // 10)]
    boxplot = [data for idx, data in enumerate(boxplot) if idx in x_values]
    marker_props = dict(marker='.', markerfacecolor='black', markersize=2, linestyle='none')
    ax1.boxplot(boxplot, widths=num_steps//20, positions=x_values, flierprops=marker_props)
    
    ax1.plot(sample_probs, linestyle='', marker='x', color='green', alpha=0.4, label=r'Best Sample $p(y \mid x)$')
    # ax1.plot(objs_exp, color='blue', alpha=0.8, label=r'Fractional $\exp \mathbb{E}[\log p(y|x)]$')
    ax1.plot(approx_prob, color='red', alpha=0.8, label=r'Sampling Approx $E[p (y  \mid x)]$')
    ax1.plot(pyx, color='orange', alpha=0.9, label=r'Integral Solution $p(y \mid x)$')

    ax1.tick_params(axis='y')
    ax1.legend(fontsize="8")

    if not os.path.exists(f"graphs/{args.folder}"):
        os.makedirs(f"graphs/{args.folder}")

    plt.title(f'id {rna_id}, init={init}, lr={learning_rate}, k={sample_size}, time={time:.2f}')
    plt.savefig(f'graphs/{args.folder}/{rna_id}.png', format="png", bbox_inches="tight")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--y", type=str, default="")
    parser.add_argument("--folder", type=str, default="temp")
    parser.add_argument("--file", type=str, default="0")
    args = parser.parse_args()

    results_path = f'./results/{args.folder}/{args.file}.txt'
    with open(results_path, 'r') as result_file:
        process_result_file(args.file, result_file, args)


if __name__ == '__main__':
    main()