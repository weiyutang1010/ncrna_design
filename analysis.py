import os, sys
import numpy as np
import argparse
import matplotlib.pyplot as plt
import concurrent.futures
from collections import defaultdict

import RNA

def prob(seq, ss, scale=True):
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

def print_subopt_result(structure, energy, data):
    ss_list = []
    if not structure == None:
        data['ss_list'].append((energy, structure))
        data['counter'] = data['counter'] + 1

def subopt(seq, e=0):
    fc = RNA.fold_compound(seq)
    fc.subopt_cb(e, print_subopt_result, subopt_data)
    subopt_data['ss_list'] = sorted(subopt_data['ss_list'])
    return subopt_data

def eval_seq(seq, ss, scale=True):
    subopt_data = { 'counter' : 0, 'sequence' : seq, 'ss_list': []}
    
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    fc.bpp()
    fc.subopt_cb(0, print_subopt_result, subopt_data)

    pr = fc.pr_structure(ss)
    ed = fc.ensemble_defect(ss)
    
    mfe_structs = [st for e, st in subopt_data['ss_list']]
    is_mfe = ss in mfe_structs
    is_umfe = is_mfe and subopt_data['counter'] == 1

    return seq, pr, ed, is_mfe, is_umfe

def graph(rna_id, lines, avg_pyx, integral_pyx, sampled_pyx, boxplot, args):
    plt.rcParams["figure.figsize"] = [7.50, 4.50]
    plt.rcParams["figure.autolayout"] = True

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Step')
    ax1.set_ylabel('Boltzmann Probability')

    n, rna_struct = len(lines[0]), lines[0]
    init = lines[1].split(', ')[1].split(': ')[1]
    learning_rate = lines[2].split(', ')[0].split(': ')[1]
    # sample_size = lines[3].split(', ')[0].split(': ')[1]
    sample_size = 2500

    time = 0.0
    if lines[-2].startswith("Total Time: "):
        time = float(lines[-2].split(': ')[1])

    # box plot
    num_steps = len(avg_pyx)
    x_values = [x for x in range(0, num_steps, num_steps // 10)]
    boxplot = [data for idx, data in enumerate(boxplot) if idx in x_values]
    marker_props = dict(marker='.', markerfacecolor='black', markersize=2, linestyle='none')
    ax1.boxplot(boxplot, widths=num_steps//20, positions=x_values, flierprops=marker_props)
    
    # ax1.plot(objs_exp, color='blue', alpha=0.8, label=r'Fractional $\exp \mathbb{E}[\log p(y|x)]$')
    ax1.plot(sampled_pyx, linestyle='', marker='x', color='green', alpha=0.4, label=r'Best Sample $p(y \mid x)$')
    ax1.plot(avg_pyx, color='red', alpha=0.8, label=r'Sampling Approx $E[p (y  \mid x)]$')
    ax1.plot(integral_pyx, color='orange', alpha=0.9, label=r'Integral Solution $p(y \mid x)$')

    ax1.tick_params(axis='y')
    ax1.legend(fontsize="8")

    if not os.path.exists(f"graphs/{args.folder}"):
        os.makedirs(f"graphs/{args.folder}")

    plt.title(f'id {rna_id}, init={init}, lr={learning_rate}, k={sample_size}, time={time:.2f}')
    plt.savefig(f'graphs/{args.folder}/{rna_id}.png', format="png", bbox_inches="tight")

def process_result_file(rna_id, result_file, args):
    lines = result_file.read().split('\n')

    if len(lines) < 5:
        exit(0)

    n, rna_struct = len(lines[0]), lines[0]

    avg_pyx = [] # avg p(y | x) of sampled sequences
    integral_seqs, integral_pyx = [], [] # integral solution at each iteration
    sampled_seqs, sampled_pyx = [], [] # best sampled solution at each iteration
    boxplot = []
    seqs = set()

    # File reading
    for idx, line in enumerate(lines):
        if line.startswith("Boxplot: "):
            values = line.split(': ')[1].split(' ')
            values = [float(value) for value in values if len(value) > 0]
            boxplot.append(values)

        if line.startswith("step:"):
            values = line.split(', ')
            seqs.add(values[3].split(': ')[1])
            avg_pyx.append(float(values[2].split(': ')[1]))
            integral_pyx.append(float(values[4].split(': ')[1]))

        if line.startswith("best samples"):
            j = idx + 1
            seqs.add(lines[j].split(' ')[0])
            sampled_pyx.append(float(lines[j].split(' ')[1]))

    if len(avg_pyx) < 1:
        exit(0)

    graph(rna_id, lines, avg_pyx, integral_pyx, sampled_pyx, boxplot, args)

    print("Seq length: ", len(seqs))
    seqs_stats = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        results = [executor.submit(eval_seq, seq, rna_struct) for seq in seqs]
        concurrent.futures.wait(results)

        for future in results:
            result = future.result()
            seqs_stats.append(result)

    print("Seq stats length: ", len(seqs_stats))

    best_pyx_solution = max(seqs_stats, key=lambda x: x[1])
    best_ned_solution = min(seqs_stats, key=lambda x: x[2])
    mfe_solutions = [stat[0] for stat in seqs_stats if stat[3]]
    umfe_solutions = [stat[0] for stat in seqs_stats if stat[4]]

    print("Best p(y|x) solution: ", best_pyx_solution[0], best_pyx_solution[1])
    print("Best NED solution: ", best_ned_solution[0], best_ned_solution[2])
    print("MFE solution: ", *mfe_solutions[:1])
    print("UMFE solution: ", *umfe_solutions[:1])

    print(f"id {rna_id} done", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--y", type=str, default="")
    parser.add_argument("--folder", type=str, default="temp")
    parser.add_argument("--file", type=str, default="0")
    args = parser.parse_args()

    results_path = f'./results/{args.folder}/{args.file}'

    # remove file extension
    rna_id = os.path.splitext(args.file)[0]

    with open(results_path, 'r') as result_file:
        process_result_file(rna_id, result_file, args)

if __name__ == '__main__':
    main()