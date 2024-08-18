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
    avg_obj = avg_obj[:2800]
    avg_pyx = avg_pyx[:2800]
    integral_pyx = integral_pyx[:2800]
    sampled_pyx = sampled_pyx[:2800]

    plt.rcParams["figure.figsize"] = [7.50, 4.50]
    # plt.rcParams["figure.figsize"] = [13, 5.50]
    plt.rcParams["figure.autolayout"] = True

    # fig, (ax1, ax2) = plt.subplots(2, 1,  gridspec_kw={'height_ratios': [2, 1]}, sharex=True)
    fig, ax1 = plt.subplots()

    # ax2.set_xlabel('Step', fontsize=fontsize)
    ax1.set_xlabel('Step', fontsize=fontsize)
    ax1.set_ylabel('Boltzmann Probability', fontsize=fontsize)

    n, rna_struct = len(lines[0]), lines[0]
    init = lines[1].split(', ')[1].split(': ')[1]
    learning_rate = lines[2].split(', ')[0].split(': ')[1]
    # sample_size = lines[3].split(', ')[0].split(': ')[1]
    sample_size = 2500

    # lr change
    if len(lr_idx) > 0:
        plt.axvline(x=lr_idx[0], color='black', linestyle='--', alpha=0.25, label='lr Decay')
        for idx in lr_idx[1:]:
            plt.axvline(x=idx, color='black', linestyle='--', alpha=0.25)
    

    # box plot
    num_steps = len(avg_pyx)
    # num_steps = 24
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

    # integral_pyx[0] = 0.020
    # boxplot = np.transpose(boxplot)

    # # entropy
    # x = [i for i in range(100)]
    # entropy = [9.5393, 9.0479, 8.9508, 8.994, 8.7523, 8.8146, 8.4899, 8.5017, 8.6286, 8.5619, 8.4882, 8.3358, 8.0615, 8.1179, 8.1122, 8.0453, 8.1063, 7.7116, 7.8385, 7.8599, 7.8379, 7.7168, 7.8763, 7.457, 7.5019, 7.5649, 7.5443, 7.4862, 7.5239, 7.5814, 7.5903, 7.4797, 7.4552, 7.4828, 7.5537, 7.6256, 7.6498, 7.4166, 7.3056, 7.4072, 7.4749, 7.2726, 7.3185, 7.3055, 7.255, 7.0351, 7.0665, 7.1164, 7.0683, 7.1253, 7.0909, 7.155, 7.1365, 7.1556, 7.1772, 7.2508, 7.1408, 6.8837, 6.9262, 7.0275, 6.9734, 6.9092, 6.9517, 6.938, 6.9857, 7.0115, 7.0723, 6.8329, 6.9316, 6.9986, 7.0109, 6.8077, 6.8454, 6.8592, 6.9708, 7.0391, 6.9526, 6.9395, 7.0049, 7.029, 6.9795, 7.069, 7.0427, 7.0925, 6.9257, 7.0033, 7.0466, 7.0361, 6.9614, 6.9534, 7.02, 6.9431, 6.9433, 7.0258, 6.7789, 6.7654, 6.8757, 6.9372, 6.904, 7.074]
    # ax2.bar(x[:21], entropy[:21])
    # ax2.set_ylim(7.5, 9.6)
    # ax2.invert_yaxis()
    # ax2.set_ylabel("Entropy", fontsize=fontsize)
    # ax2.set_yticks([7.5, 8.0, 8.5, 9.0, 9.5])

    # ax1.plot(sampled_pyx[:21], linestyle='', marker='o', markerfacecolor='None', markeredgecolor='green',  alpha=0.5, ms=11)
    # for idx, data in enumerate(boxplot):
    #     if idx == 0:
    #         ax1.plot(data[:21], linestyle='', marker='x', color='green', alpha=0.3, label=r'samples (best in circle)')
    #     else:
    #         ax1.plot(data[:21], linestyle='', marker='x', color='green', alpha=0.3)

    # ax1.plot(integral_pyx[:21], color='orange', alpha=0.9, label=r'integral solution')
    # ax1.plot(avg_pyx[:21], color='red', alpha=0.8, label=r'arith. mean')
    # ax1.plot(np.exp(-1 * np.array(avg_obj[:21])), color='blue', alpha=0.8, label=r'geom. mean $(e^{-\mathcal{J}})$')

    # ax1.set_xticks([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    # ax1.set_xlim(-0.5, 20.5)

    # ax1.set_xlim(0, 3000)
    ax1.set_xlim(-100, 2900)
    ax1.axvline(x=2000, color='black', linestyle='--', alpha=0.5)
    ax1.tick_params(labelbottom=True, axis='both', which='major', labelsize=13)
    # ax2.tick_params(labelbottom=False, axis='both', which='major', labelsize=13)
    # plt.subplots_adjust(wspace=0, hspace=0)

    ax1.legend(fontsize=fontsize)
    # ax1.legend(fontsize=fontsize, loc='upper right', bbox_to_anchor=(1.36, 1.04))

    if not os.path.exists(f"graphs/{args.folder}"):
        os.makedirs(f"graphs/{args.folder}")

    save_path = f'graphs/{args.folder}/{rna_id}.pdf'
    # plt.title(f'id {rna_id}, init={init}, lr={learning_rate}, k={sample_size}, time={time:.2f}')
    # plt.title(f'Learning Curves, (((...)))')
    plt.savefig(save_path, bbox_inches="tight")
    # plt.savefig(save_path, format="png", bbox_inches="tight")
    print(f"Puzzle {rna_id} saved to {save_path}", file=sys.stderr)

def graph_ned(rna_id, lines, avg_ned, integral_ned, sampled_ned, boxplot, lr_idx, args):
    plt.rcParams["figure.figsize"] = [7.50, 4.50]
    plt.rcParams["figure.autolayout"] = True

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Step')
    ax1.set_ylabel('log NED')

    n, rna_struct = len(lines[0]), lines[0]
    init = lines[1].split(', ')[1].split(': ')[1]
    learning_rate = lines[2].split(', ')[0].split(': ')[1]
    # sample_size = lines[3].split(', ')[0].split(': ')[1]
    sample_size = 2500

    time = 0.0
    if lines[-2].startswith("Total Time: "):
        time = float(lines[-2].split(': ')[1])

    # lr change
    if len(lr_idx) > 0:
        plt.axvline(x=lr_idx[0], color='black', linestyle='--', alpha=0.25, label='lr Decay')
        for idx in lr_idx[1:]:
            plt.axvline(x=idx, color='black', linestyle='--', alpha=0.25)

    # box plot
    num_steps = len(avg_ned)
    x_values = [x for x in range(0, num_steps, num_steps // 10)]
    boxplot = [data for idx, data in enumerate(boxplot) if idx in x_values]
    marker_props = dict(marker='.', markerfacecolor='black', markersize=2, linestyle='none')
    ax1.boxplot(boxplot, widths=num_steps//20, positions=x_values, flierprops=marker_props)
    
    # ax1.plot(objs_exp, color='blue', alpha=0.8, label=r'Fractional $\exp \mathbb{E}[\log p(y|x)]$')
    ax1.plot(sampled_ned, linestyle='', marker='x', color='green', alpha=0.4, label=r'Best Sample NED')
    ax1.plot(avg_ned, color='red', alpha=0.8, label=r'Average NED')
    ax1.plot(integral_ned, color='orange', alpha=0.9, label=r'Integral Solution NED')

    ax1.tick_params(axis='y')
    ax1.legend(fontsize="8")

    if not os.path.exists(f"graphs/{args.folder}"):
        os.makedirs(f"graphs/{args.folder}")

    save_path = f'graphs/{args.folder}/{rna_id}.png'
    plt.title(f'id {rna_id}, init={init}, lr={learning_rate}, k={sample_size}, time={time:.2f}')
    plt.savefig(save_path, format="png", bbox_inches="tight")
    print(f"Puzzle {rna_id} saved to {save_path}", file=sys.stderr)

def graph_dist(rna_id, lines, avg_dist, integral_dist, sampled_dist, boxplot, lr_idx, args):
    # plt.rcParams["text.usetex"] = True
    plt.rcParams["figure.figsize"] = [7.50, 4.50]
    plt.rcParams["figure.autolayout"] = True

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Step')
    ax1.set_ylabel('- Nemo Objective')

    n, rna_struct = len(lines[0]), lines[0]
    init = lines[1].split(', ')[1].split(': ')[1]
    learning_rate = lines[2].split(', ')[0].split(': ')[1]
    # sample_size = lines[3].split(', ')[0].split(': ')[1]
    sample_size = 2500

    time = 0.0
    if lines[-2].startswith("Total Time: "):
        time = float(lines[-2].split(': ')[1])

    # lr change
    if len(lr_idx) > 0:
        plt.axvline(x=lr_idx[0], color='black', linestyle='--', alpha=0.25, label='lr Decay')
        for idx in lr_idx[1:]:
            plt.axvline(x=idx, color='black', linestyle='--', alpha=0.25)

    plt.axhline(y=-1.0, color='blue', linestyle='--', alpha=0.25, label='mfe found')

    # box plot
    num_steps = len(avg_dist)
    x_values = [x for x in range(0, num_steps, num_steps // 10)]
    boxplot = [data for idx, data in enumerate(boxplot) if idx in x_values]
    marker_props = dict(marker='.', markerfacecolor='black', markersize=2, linestyle='none')
    ax1.boxplot(boxplot, widths=num_steps//20, positions=x_values, flierprops=marker_props)
    
    # ax1.plot(objs_exp, color='blue', alpha=0.8, label=r'Fractional $\exp \mathbb{E}[\log p(y|x)]$')
    ax1.plot(sampled_dist, linestyle='', marker='x', color='green', alpha=0.4, label=r'Best Sample')
    ax1.plot(avg_dist, color='red', alpha=0.8, label=r'Average')
    ax1.plot(integral_dist, color='orange', alpha=0.9, label=r'Integral Solution')

    ax1.tick_params(axis='y')
    ax1.legend(fontsize="8")

    if not os.path.exists(f"graphs/{args.folder}"):
        os.makedirs(f"graphs/{args.folder}")

    save_path = f'graphs/{args.folder}/{rna_id}.png'
    plt.title(f'id {rna_id}, init={init}, lr={learning_rate}, k={sample_size}, time={time:.2f}')
    plt.savefig(save_path, format="png", bbox_inches="tight")
    print(f"Puzzle {rna_id} saved to {save_path}", file=sys.stderr)

def process_result_file(rna_id, result_file, args):
    lines = result_file.read().split('\n')

    if len(lines) < 5:
        exit(0)

    n, rna_struct = len(lines[0]), lines[0]
    is_ned = 'ned' in lines[1].split(', ')[0].split(': ')[1]
    is_dist = 'dist' in lines[1].split(', ')[0].split(': ')[1]

    obj, avg_pyx = [], [] # avg p(y | x) of sampled sequences
    integral_seqs, integral_obj = [], [] # integral solution at each iteration
    sampled_seqs, sampled_obj = [], [] # best sampled solution at each iteration
    boxplot = []
    prev_lr = float(lines[2].split(', ')[0].split(': ')[1])
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
            seq = values[3].split(': ')[1]

            seqs.add(seq)
            obj.append(float(values[1].split(': ')[1]))
            avg_pyx.append(float(values[2].split(': ')[1]))
            integral_obj.append(float(values[4].split(': ')[1]))

            if seq not in seq_step:
                seq_step[seq] = step

            if values[5].split(': ')[0] == 'learning rate':
                lr = float(values[5].split(': ')[1])
                if lr != prev_lr:
                    prev_lr = lr
                    lr_idx.append(step)

        if line.startswith("best samples"):
            j = idx + 1
            seq = lines[j].split(' ')[0]
            seq_obj = float(lines[j].split(' ')[1])

            seqs.add(seq)
            sampled_obj.append(seq_obj)
            if seq not in seq_step:
                seq_step[seq] = step


    if len(avg_pyx) < 1:
        exit(0)

    if is_ned:
        graph_ned(rna_id, lines, obj, integral_obj, sampled_obj, boxplot, lr_idx, args)
    elif is_dist:
        graph_dist(rna_id, lines, obj, integral_obj, sampled_obj, boxplot, lr_idx, args)
    else:
        graph_prob(rna_id, lines, obj, avg_pyx, integral_obj, sampled_obj, boxplot, lr_idx, args)

    exit(0)

    seqs_stats = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
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
    best_ediff_solution = min(seqs_stats, key=lambda x: x[6])

    print("Best p(y|x) solution: ", best_pyx_solution[0], best_pyx_solution[1], seq_step[best_pyx_solution[0]])
    print("Best NED solution: ", best_ned_solution[0], best_ned_solution[2], seq_step[best_ned_solution[0]])
    print("MFE solution: ", len(mfe_solutions), *mfe_solutions[:1])
    print("UMFE solution: ", len(umfe_solutions), *umfe_solutions[:1])
    print("Best d(MFE(x), y) solution: ", best_dist_solution[0], best_dist_solution[5], seq_step[best_dist_solution[0]])
    print("Best free energy diff solution: ", best_ediff_solution[0], best_ediff_solution[6], seq_step[best_ediff_solution[0]])

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