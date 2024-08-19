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

def e_diff(seq, ss):
    ss_mfe = mfe(seq)[0]
    return abs(energy(seq, ss_mfe) - energy(seq, ss))

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

    dist = structural_dist(seq, ss)
    energy_diff = e_diff(seq, ss)

    return seq, pr, ed, is_mfe, is_umfe, dist, energy_diff

def graph_prob(rna_id, lines, avg_obj, avg_pyx, integral_pyx, sampled_pyx, boxplot, lr_idx, args, fontsize=14):
    plt.rcParams["figure.figsize"] = [9.50, 4.50]
    plt.rcParams["figure.autolayout"] = True

    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Step', fontsize=fontsize)
    ax1.set_ylabel('Boltzmann Probability', fontsize=fontsize)

    n, rna_struct = len(lines[0]), lines[0]
    init = lines[1].split(', ')[1].split(': ')[1]
    learning_rate = lines[3].split(', ')[0].split(': ')[1]
    sample_size = lines[6].split(', ')[0].split(': ')[1]

    # lr change
    if len(lr_idx) > 0:
        plt.axvline(x=lr_idx[0], color='black', linestyle='--', alpha=0.25, label='lr decay')
        for idx in lr_idx[1:]:
            plt.axvline(x=idx, color='black', linestyle='--', alpha=0.25)
    
    # box plot
    if len(boxplot) > 0:
        num_steps = len(avg_pyx)
        x_values = [x for x in range(0, num_steps, (num_steps + 9) // 10)]
        boxplot = [data for idx, data in enumerate(boxplot) if idx in x_values]
        marker_props = dict(marker='.', markerfacecolor='black', markersize=2, linestyle='none')
        ax1.boxplot(boxplot, widths=num_steps//20, positions=x_values, flierprops=marker_props)

    # calculate geom. mean
    objs_exp = np.exp(-1 * np.array(avg_obj))
    
    ax1.plot(sampled_pyx, linestyle='', marker='o', markerfacecolor='None', color='green', alpha=0.3, label=r'best sample')
    ax1.plot(integral_pyx, color='orange', alpha=0.9, label=r'integral solution')
    ax1.plot(avg_pyx, color='red', alpha=0.8, label=r'arith. mean')
    ax1.plot(objs_exp, color='blue', alpha=0.8, label=r'geom. mean $(e^{-\mathcal{J}})$')

    ax1.tick_params(labelbottom=True, axis='both', which='major', labelsize=13)

    ax1.legend(fontsize=fontsize)

    if not os.path.exists(f"graphs/{args.folder}"):
        os.makedirs(f"graphs/{args.folder}")

    save_path = f'graphs/{args.folder}/{rna_id}.pdf'
    # plt.title(f'id {rna_id}, init={init}, sample_size={sample_size}')
    plt.savefig(save_path, bbox_inches="tight")
    print(f"Puzzle {rna_id} saved to {save_path}", file=sys.stderr)

def graph(rna_id, objective, lines, avg_obj, integral, sampled, boxplot, lr_idx, args, fontsize=14):
    plt.rcParams["figure.figsize"] = [9.50, 4.50]
    plt.rcParams["figure.autolayout"] = True

    fig, ax1 = plt.subplots()

    ax1.set_xlabel('Step', fontsize=fontsize)
    if objective == 'ned':
        ax1.set_ylabel('Normalized Ensemble Defect', fontsize=fontsize)
    elif objective == 'dist':
        ax1.set_ylabel('Structural Distance', fontsize=fontsize)
    elif objective == 'ddg':
        ax1.set_ylabel('Free Energy Gap', fontsize=fontsize)
    else:
        ax1.set_ylabel(objective, fontsize=fontsize)

    n, rna_struct = len(lines[0]), lines[0]
    init = lines[1].split(', ')[1].split(': ')[1]
    learning_rate = lines[3].split(', ')[0].split(': ')[1]
    sample_size = lines[6].split(', ')[0].split(': ')[1]

    # lr change
    if len(lr_idx) > 0:
        plt.axvline(x=lr_idx[0], color='black', linestyle='--', alpha=0.25, label='lr decay')
        for idx in lr_idx[1:]:
            plt.axvline(x=idx, color='black', linestyle='--', alpha=0.25)
    
    # box plot
    if len(boxplot) > 0:
        num_steps = len(avg_obj)
        x_values = [x for x in range(0, num_steps, (num_steps + 9) // 10)]
        boxplot = [data for idx, data in enumerate(boxplot) if idx in x_values]
        marker_props = dict(marker='.', markerfacecolor='black', markersize=2, linestyle='none')
        ax1.boxplot(boxplot, widths=num_steps//20, positions=x_values, flierprops=marker_props)

    # calculate geom. mean
    
    ax1.plot(sampled, linestyle='', marker='o', markerfacecolor='None', color='green', alpha=0.3, label=r'best sample')
    ax1.plot(integral, color='orange', alpha=0.9, label=r'integral solution')
    ax1.plot(avg_obj, color='blue', alpha=0.8, label=r'geom. mean $(e^{-\mathcal{J}})$')

    ax1.tick_params(labelbottom=True, axis='both', which='major', labelsize=13)
    ax1.legend(fontsize=fontsize)

    if not os.path.exists(f"graphs/{args.folder}"):
        os.makedirs(f"graphs/{args.folder}")

    save_path = f'graphs/{args.folder}/{rna_id}.pdf'
    # plt.title(f'id {rna_id}, init={init}, sample_size={sample_size}')
    plt.savefig(save_path, bbox_inches="tight")
    print(f"Puzzle {rna_id} saved to {save_path}", file=sys.stderr)

def process_result_file(rna_id, result_file, args):
    lines = result_file.read().split('\n')

    if len(lines) < 10:
        exit(0)

    n, rna_struct = len(lines[0]), lines[0]
    objective = lines[1].split(', ')[0].split(': ')[1]

    obj, avg_pyx = [], [] # avg p(y | x) of sampled sequences
    integral_seqs, integral_obj = [], [] # integral solution at each iteration
    sampled_seqs, sampled_obj = [], [] # best sampled solution at each iteration
    boxplot = []
    prev_lr = float(lines[3].split(', ')[0].split(': ')[1])
    lr_idx = [] # track when does lr changes
    seqs = set()
    seq_step = {}

    # File reading
    for idx, line in enumerate(lines):
        if line.startswith("Boxplot: "):
            values = line.split(': ')[1].split(' ')
            try:
                values = [float(value) for value in values if len(value) > 0]
            except ValueError as e:
                print(f"Puzzle {rna_id} Error: {e}", file=sys.stderr)
                continue

            boxplot.append(values)

        if line.startswith("step:"):
            values = line.split(', ')
            step = int(values[0].split(': ')[1])

            obj.append(float(values[1].split(': ')[1]))
            avg_pyx.append(float(values[2].split(': ')[1]))

            lr = float(values[3].split(': ')[1])
            if lr != prev_lr:
                prev_lr = lr
                lr_idx.append(step)

        if line.startswith("integral solution"):
            seq = line.split(': ')[1].split(' ')[0]
            seq_obj = float(line.split(': ')[1].split(' ')[1])

            seqs.add(seq)
            integral_obj.append(seq_obj)
            if seq not in seq_step:
                seq_step[seq] = step

        if line.startswith("best samples"):
            j = idx + 1 # evaluate the best sample
            seq = lines[j].split(' ')[0]
            seq_obj = float(lines[j].split(' ')[1])

            seqs.add(seq)
            sampled_obj.append(seq_obj)
            if seq not in seq_step:
                seq_step[seq] = step


    if len(avg_pyx) < 1:
        exit(0)

    # graph solutions
    if not args.no_graph:
        if objective == "prob":
            graph_prob(rna_id, lines, obj, avg_pyx, integral_obj, sampled_obj, boxplot, lr_idx, args)
        else:
            graph(rna_id, objective, lines, obj, integral_obj, sampled_obj, boxplot, lr_idx, args)

    # use vienna to reevaluate all sequences
    if not args.no_eval_seq:
        seqs_stats = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            results = [executor.submit(eval_seq, seq, rna_struct) for seq in seqs]
            concurrent.futures.wait(results)

            for future in results:
                result = future.result()
                seqs_stats.append(result)

        print("Total steps: ", len(obj))
        print("Number of unique sequence: ", len(seqs_stats))

        best_pyx_solution = max(seqs_stats, key=lambda x: x[1])
        best_ned_solution = min(seqs_stats, key=lambda x: x[2])
        mfe_solutions = [stat[0] for stat in seqs_stats if stat[3]]
        umfe_solutions = [stat[0] for stat in seqs_stats if stat[4]]
        best_dist_solution = min(seqs_stats, key=lambda x: x[5])
        best_ddg_solution = min(seqs_stats, key=lambda x: x[6])

        mfe_solution = "" if len(mfe_solutions) == 0 else mfe_solutions[0]
        umfe_solution = "" if len(umfe_solutions) == 0 else umfe_solutions[0]
        mfe_step = "" if len(mfe_solutions) == 0 else seq_step[mfe_solution]
        umfe_step = "" if len(umfe_solutions) == 0 else seq_step[umfe_solution]

        print("Format: <sequence> <value> <earliest step found>")
        print("Best p(y*|x) solution: ", best_pyx_solution[0], best_pyx_solution[1], seq_step[best_pyx_solution[0]])
        print("Best NED(x, y*) solution: ", best_ned_solution[0], best_ned_solution[2], seq_step[best_ned_solution[0]])
        print("Best d(MFE(x), y) solution: ", best_dist_solution[0], best_dist_solution[5], seq_step[best_dist_solution[0]])
        print("Best DeltaDeltaG(x, y*) solution: ", best_ddg_solution[0], best_ddg_solution[6], seq_step[best_ddg_solution[0]])
        print("MFE/uMFE Format: <number of mfe solutions> [solution 1] [earliest step found]")
        print("MFE solution: ", len(mfe_solutions), mfe_solution, mfe_step) # only print 1 mfe solution as example
        print("UMFE solution: ", len(umfe_solutions), umfe_solution, umfe_step)

    print(f"id {rna_id} done", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--y", type=str, default="")
    parser.add_argument("--folder", type=str)
    parser.add_argument("--file", type=str)
    parser.add_argument("--max_workers", type=int, default=None)
    parser.add_argument("--no_eval_seq", action="store_true", default=False, help='draw graph only')
    parser.add_argument("--no_graph", action="store_true", default=False, help='eval seqs only')
    args = parser.parse_args()

    results_path = f'./results/{args.folder}/{args.file}'

    # remove file extension
    rna_id = os.path.splitext(args.file)[0]

    with open(results_path, 'r') as result_file:
        process_result_file(rna_id, result_file, args)

if __name__ == '__main__':
    main()