import os, sys
import numpy as np
from scipy.stats import gmean
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

import concurrent.futures
import tikzplotlib
import RNA

ned = False
if len(sys.argv) > 1:
    ned = (sys.argv[1] == 'ned')
    dist = (sys.argv[1] == 'dist')
    ddg = (sys.argv[1] == 'ddg')

ids = []
puzzles = []

with open('../data/eterna/eterna100.txt', 'r') as data:
    lines = data.read().split('\n')
    ids = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
    puzzles = [line.split(' ')[1] for line in lines if len(line) > 0]

undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]
unknown_ids = [68, 97, 100]

def get_boltz_prob(seq, ss, scale=True):
    """viennaRNA boltzmann probability"""
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    Q = fc.pf()
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

# for file in ["temp.txt"]:
#     with open(file, 'r') as f:
#         lines = f.read().split('\n')

#     with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
#         futures = [executor.submit(structural_dist, line, puzzles[idx]) for idx, line in enumerate(lines)]

#     concurrent.futures.wait(futures)
#     obj_arr = [future.result() for future in futures]
#     print(np.mean(obj_arr))
# exit(0)


def parse(file_path, ned=False, dist=False, ddg=False):
    res, res1, res2, res3 = [], [], [], []
    with open(file_path, 'r') as f:
        lines = f.read().split('\n')

        obj_arr = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=64) as executor:
            if ned:
                futures = [executor.submit(ensemble_defect, line, puzzles[idx]) for idx, line in enumerate(lines)]
            elif ddg:
                futures = [executor.submit(e_diff, line, puzzles[idx]) for idx, line in enumerate(lines)]
            else:
                futures = [executor.submit(get_boltz_prob, line, puzzles[idx]) for idx, line in enumerate(lines)]
            
            concurrent.futures.wait(futures)
            obj_arr = [future.result() for future in futures]

        idx = 0
        for p_id, line in zip(ids, lines):
            obj = obj_arr[idx]

            if p_id in undesignable_ids:
                res2.append(obj)
            elif p_id in unknown_ids:
                res3.append(obj)
            else:
                res1.append(obj)

            res.append(obj)
            idx += 1

    return res, res1, res2, res3

if __name__ == '__main__':
    if ned:
        sampling_soft, sampling_soft_1, sampling_soft_2, sampling_soft_3 = parse("./softmax_new_ned.txt", ned=True)
        samfeo_5k, samfeo_5k_1,samfeo_5k_2, samfeo_5k_3 = parse("./samfeo_5k_ned.txt", ned=True)
        # samfeo_10k, samfeo_10k_1,samfeo_10k_2, samfeo_10k_3 = parse("./samfeo_10k_ned.txt", ned=True)
    elif ddg:
        sampling_soft, sampling_soft_1, sampling_soft_2, sampling_soft_3 = parse("./softmax_new_ddg.txt", ddg=True)
        samfeo_5k, samfeo_5k_1,samfeo_5k_2, samfeo_5k_3 = parse("./samfeo_5k_ddg.txt", ddg=True)
    else:
        sampling_soft, sampling_soft_1, sampling_soft_2, sampling_soft_3 = parse("./softmax_new.txt")
        samfeo_5k, samfeo_5k_1,samfeo_5k_2, samfeo_5k_3 = parse("./samfeo_5k.txt")
        # samfeo_10k, samfeo_10k_1,samfeo_10k_2, samfeo_10k_3 = parse("./samfeo_10k.txt")

    a, b, c, d = 12, 400, 1.5, 12
    puzzle_lens_1 = [len(puzzle) for p_id, puzzle in zip(ids, puzzles) if p_id not in undesignable_ids and p_id not in unknown_ids]
    map_len_1 = [c + (((x - a) * (d - c)) / (b - a)) for x in puzzle_lens_1]
    puzzle_lens_2 = [len(puzzle) for p_id, puzzle in zip(ids, puzzles) if p_id in undesignable_ids]
    map_len_2 = [c + (((x - a) * (d - c)) / (b - a)) for x in puzzle_lens_2]
    puzzle_lens_3 = [len(puzzle) for p_id, puzzle in zip(ids, puzzles) if p_id in unknown_ids]
    map_len_3 = [c + (((x - a) * (d - c)) / (b - a)) for x in puzzle_lens_3]

    id_1 = [p_id for p_id, puzzle in zip(ids, puzzles) if p_id not in undesignable_ids and p_id not in unknown_ids]
    id_2 = [p_id for p_id, puzzle in zip(ids, puzzles) if p_id in undesignable_ids]
    id_3 = [p_id for p_id, puzzle in zip(ids, puzzles) if p_id in unknown_ids]

    arr = [(12, 36), (38, 57), (61, 75), (80, 98), (100, 104), (105, 111), (116, 192), (200, 316), (337, 387), (389, 400)]
    samp = [[] for _ in range(10)]
    samf = [[] for _ in range(10)]
    for idx, x in enumerate(sampling_soft):
        # if ids[idx] in undesignable_ids:
        #     continue

        for i, bound in enumerate(arr):
            if len(puzzles[idx]) >= bound[0] and len(puzzles[idx]) <= bound[1]:
                samp[i].append(x)
                samf[i].append(samfeo_5k[idx])

    for i in range(10):
        print(f"({arr[i][0]}--{arr[i][1]},{np.mean(samp[i]):.4f})", end=" ")
    print()

    for i in range(10):
        print(f"({arr[i][0]}--{arr[i][1]},{np.mean(samf[i]):.4f})", end=" ")
    print()

    # for pid, x, y, datasize in zip(id_1, samfeo_5k_1, sampling_soft_1, map_len_1):
    #     # print(x, 4), y, 4), datasize, 6), "Designable")
    #     print(f"{pid} {x:.4f} {y:.4f} {datasize:.6f} Designable")

    # for pid, x, y, datasize in zip(id_2, samfeo_5k_2, sampling_soft_2, map_len_2):
    #     # print(x, 4), y, 4), datasize, 6), "Undesignable")
    #     print(f"{pid} {x:.4f} {y:.4f} {datasize:.6f} Undesignable")

    # for pid, x, y, datasize in zip(id_3, samfeo_5k_3, sampling_soft_3, map_len_3):
    #     # print(x, 4), y, 4), round(datasize, 6), "Unknown")
    #     print(f"{pid} {x:.4f} {y:.4f} {datasize:.6f} Unknown")



    # sampling_soft = [x for x, puzzle, p_id in zip(sampling_soft, puzzles, ids) if len(puzzle) <= 104]
    # samfeo_5k = [x for x, puzzle, p_id in zip(samfeo_5k, puzzles, ids) if len(puzzle) <= 104]
    # samfeo_10k = [x for x, puzzle, p_id in zip(samfeo_10k, puzzles, ids) if len(puzzle) <= 104]

    # print("Arithmetic Mean")
    # print("Sampling (softmax)      : ", np.mean(sampling_soft))
    # print("SAMFEO (5k)             : ", np.mean(samfeo_5k))
    # # print("SAMFEO (10k)            : ", np.mean(samfeo_10k))
    # print()


    exit(0)
    # remove puzzle 22 to calculate geometric mean
    if ned:
        idx = sampling_soft.index(0.0)
        del sampling_soft[idx]
        del samfeo_5k[idx]
        del samfeo_10k[idx]

        del ids[idx]
        del puzzles[idx]

        idx = sampling_soft_1.index(0.0)
        del sampling_soft_1[idx]
        del samfeo_5k_1[idx]
        del samfeo_10k_1[idx]

    print("Geometric Mean")
    print("Sampling (softmax)      : ", gmean(sampling_soft))
    print("SAMFEO (5k)             : ", gmean(samfeo_5k))
    print("SAMFEO (10k)            : ", gmean(samfeo_10k))
    print()

    print("Geometric Mean (w/o undesignable)")
    print("Sampling (softmax)      : ", gmean(sampling_soft_1 + sampling_soft_3))
    print("SAMFEO (5k)             : ", gmean(samfeo_5k_1 + samfeo_5k_3))
    print("SAMFEO (10k)            : ", gmean(samfeo_10k_1 + samfeo_10k_3))
    print()

    exit(0)

    # puzzles = [a for a, b in zip(puzzles, ids) if b not in undesignable_ids]

    def bar_plot(sampling_soft, samfeo_5k, samfeo_10k, puzzles, undesignable_ids, ylabel, title, save_path, ned=False, geom=False, no_undesignable=False):
        if no_undesignable:
            puzzles = [a for a, b in zip(puzzles, ids) if b not in undesignable_ids]
            sampling_soft = [a for a, b in zip(sampling_soft, ids) if b not in undesignable_ids]
            samfeo_5k = [a for a, b in zip(samfeo_5k, ids) if b not in undesignable_ids]
            samfeo_10k = [a for a, b in zip(samfeo_10k, ids) if b not in undesignable_ids]

        if geom:
            start = 9 if no_undesignable else 9
            skip = 9 if no_undesignable else 10
        else:
            start = 10
            skip = 8 if no_undesignable else 10
            
        x = np.arange(start, len(puzzles) + 1, skip)

        sampling_bar = []
        samfeo_bar = []
        sampling_proj_bar = []

        prev = 0
        xticks_labels = []
        lengths = [len(s) for s in puzzles]

        for i in range(start, len(puzzles) + 1, skip):
            if geom:
                sampling_bar.append(gmean(sampling_soft[prev:i]))
                samfeo_bar.append(gmean(samfeo_5k[prev:i]))
                sampling_proj_bar.append(gmean(samfeo_10k[prev:i]))
            else:
                sampling_bar.append(np.mean(sampling_soft[prev:i]))
                samfeo_bar.append(np.mean(samfeo_5k[prev:i]))
                sampling_proj_bar.append(np.mean(samfeo_10k[prev:i]))

            xticks_labels.append(f'{min(lengths[prev: i])} - {max(lengths[prev: i])}')
            prev = i

        width = 2.0

        fig, ax = plt.subplots(figsize=(14, 6))
        bars1 = ax.bar(x - width, sampling_bar, width, label='Sampling (Softmax, Adam)', color='lightblue', edgecolor='black')
        bars2 = ax.bar(x, samfeo_bar, width, label='SAMFEO (5k steps)', color='lightgreen', edgecolor='black')
        bars3 = ax.bar(x + width, sampling_proj_bar, width, label='SAMFEO (10k steps)', color='lightcoral', edgecolor='black')

        # xticks_labels = [f'{i-9} - {i}' for i in x]

        # Add labels and title
        # ax.set_ylabel('Average p(y | x)')
        # ax.set_title('Average p(y | x) grouped by length')
        # ax.set_ylabel('Geometric Average p(y | x)')
        # ax.set_title('Geometric p(y | x) grouped by length')

        ax.set_xlabel('Puzzles Length Range')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(xticks_labels)
        plt.legend()

        plt.savefig(save_path, dpi=400, bbox_inches='tight')
        plt.close()
        print(f"Figure saved to {save_path}")

    def scatter_plot(data1, data1_2, data1_3, data2, data2_2, data2_3, marker_size_1, marker_size_2, marker_size_3, xlabel, ylabel, title, save_file, use_log=False, ned=False):
        if ned:
            data1 = 1 - np.array(data1)
            data1_2 = 1 - np.array(data1_2)
            data1_3 = 1 - np.array(data1_3)
            data2 = 1 - np.array(data2)
            data2_2 = 1 - np.array(data2_2)
            data2_3 = 1 - np.array(data2_3)

        # Create a scatter plot
        if use_log:
            if ned:
                fig, ax = plt.subplots()
            else:
                fig, ax = plt.subplots(figsize=(18,6))

            ax.scatter(np.log(data1), np.log(data2), s=marker_size_1, c='blue', alpha=0.4, marker='x', label='Designable')
            ax.scatter(np.log(data1_2), np.log(data2_2), s=marker_size_2, c='red', alpha=0.4, marker='x', label='Undesignable')
            ax.scatter(np.log(data1_3), np.log(data2_3), s=marker_size_3, c='green', alpha=0.4, marker='x', label='Unknown')

            if ned:
                a, b = -0.45, 0.05
                ax.plot([a, b], [a, b], color='black', linestyle='--', alpha=0.5, linewidth=1)
                ax.set_ylim(a, b)
                ax.set_xlim(a, b)
                pass
            else:
                min_point = min(min(np.log(data1_2)), min(np.log(data2_2)))
                ax.plot([min_point, 3], [min_point, 3], color='black', linestyle='--', alpha=0.5, linewidth=1)
                ax.set_ylim(-45, 3)
                ax.set_xlim(-125, 3)
        else:
            fig, ax = plt.subplots()
            ax.scatter(data1, data2, s=marker_size_1, c='blue', alpha=0.4, marker='x', label='Designable')
            ax.scatter(data1_2, data2_2, s=marker_size_2, c='red', alpha=0.4, marker='x', label='Undesignable')
            ax.scatter(data1_3, data2_3, s=marker_size_3, c='green', alpha=0.4, marker='x', label='Unknown')

            if ned:
                ax.plot([0.6, 1.05], [0.6, 1.05], color='black', linestyle='--', alpha=0.5, linewidth=1)
            else:
                ax.plot([-0.05, 1], [-0.05, 1], color='black', linestyle='--', alpha=0.5, linewidth=1)

            if ned:
                ax.set_xlim(0.6, 1.05)
                ax.set_ylim(0.6, 1.05)
            else:
                ax.set_xlim(-0.05, 1.0)
                ax.set_ylim(-0.05, 1.0)


        # Create a scatter plot
        if not use_log and not ned:
            axins = zoomed_inset_axes(ax, 1.5, loc=4) # zoom = 6
            axins.scatter(data1, data2, s=marker_size_1, c='blue', alpha=0.4, marker='x', label='Designable')
            axins.scatter(data1_2, data2_2, s=marker_size_2, c='red', alpha=0.4, marker='x', label='Undesignable')
            axins.scatter(data1_3, data2_3, s=marker_size_3, c='green', alpha=0.4, marker='x', label='Unknown')

            axins.plot([-0.05, 1], [-0.05, 1], color='black', linestyle='--', alpha=0.5, linewidth=1)

            # sub region of the original image
            x1, x2, y1, y2 = -0.05, 0.2, -0.05, 0.2
            axins.set_xlim(x1, x2)
            axins.set_ylim(y1, y2)

            plt.xticks(visible=False)
            plt.yticks(visible=False)
            plt.tick_params(left = False, bottom = False) 

            # draw a bbox of the region of the inset axes in the parent axes and
            # connecting lines between the bbox and the inset axes area
            mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

        # Add labels and title
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc=2)

        plt.savefig(f"{save_file}.png", dpi=400, bbox_inches='tight')
        print(f"Scatter plot saved to {save_file}.png")
        tikzplotlib.save(f"{save_file}.tex")
        print(f"Scatter plot saved to {save_file}.tex")
        plt.close()


    # a, b, c, d = 12, 400, 15, 150
    # puzzle_lens_1 = [len(puzzle) for p_id, puzzle in zip(ids, puzzles) if p_id not in undesignable_ids and p_id not in unknown_ids]
    # map_len_1 = [c + (((x - a) * (d - c)) / (b - a)) for x in puzzle_lens_1]
    # puzzle_lens_2 = [len(puzzle) for p_id, puzzle in zip(ids, puzzles) if p_id in undesignable_ids]
    # map_len_2 = [c + (((x - a) * (d - c)) / (b - a)) for x in puzzle_lens_2]
    # puzzle_lens_3 = [len(puzzle) for p_id, puzzle in zip(ids, puzzles) if p_id in unknown_ids]
    # map_len_3 = [c + (((x - a) * (d - c)) / (b - a)) for x in puzzle_lens_3]

    if not ned:
        # scatter_plot(samfeo_5k_1, samfeo_5k_2, samfeo_5k_3, sampling_soft_1, sampling_soft_2, sampling_soft_3, puzzle_lens_1, puzzle_lens_2, puzzle_lens_3, "SAMFEO", "Sampling (Softmax)", "p(y | x) of Sampling vs. SAMFEO (5k)", "samfeo_5k")
        # scatter_plot(samfeo_10k_1, samfeo_10k_2, samfeo_10k_3, sampling_soft_1, sampling_soft_2, sampling_soft_3, puzzle_lens_1, puzzle_lens_2, puzzle_lens_3, "SAMFEO (10k steps)", "Sampling (softmax, adam)", "p(y | x) of Sampling vs. SAMFEO (10k)", "samfeo_10k")
        # scatter_plot(samfeo_5k_1, samfeo_5k_2, samfeo_5k_3, sampling_soft_1, sampling_soft_2, sampling_soft_3, puzzle_lens_1, puzzle_lens_2, puzzle_lens_3, "SAMFEO", "Sampling (Softmax)", "log p(y | x) of Sampling vs. SAMFEO (5k)", "samfeo_5k_log", use_log=True)
        # scatter_plot(samfeo_10k_1, samfeo_10k_2, samfeo_10k_3, sampling_soft_1, sampling_soft_2, sampling_soft_3, puzzle_lens_1, puzzle_lens_2, puzzle_lens_3, "SAMFEO (10k steps)", "Sampling (softmax, adam)", "log p(y | x) of Sampling vs. SAMFEO (10k)", "samfeo_10k_log", use_log=True)
        pass

    # NED plot
    if ned:
        # scatter_plot(samfeo_5k_1, samfeo_5k_2, samfeo_5k_3, sampling_soft_1, sampling_soft_2, sampling_soft_3, puzzle_lens_1, puzzle_lens_2, puzzle_lens_3, "SAMFEO (5k steps)", "Sampling (softmax, adam)", "1 - NED of Sampling vs. SAMFEO (5k)", "samfeo_5k_ned", ned=True)
        # scatter_plot(samfeo_10k_1, samfeo_10k_2, samfeo_10k_3, sampling_soft_1, sampling_soft_2, sampling_soft_3, puzzle_lens_1, puzzle_lens_2, puzzle_lens_3, "SAMFEO (10k steps)", "Sampling (softmax, adam)", "1 - NED of Sampling vs. SAMFEO (10k)", "samfeo_10k_ned", ned=True)
        # scatter_plot(samfeo_5k_1, samfeo_5k_2, samfeo_5k_3, sampling_soft_1, sampling_soft_2, sampling_soft_3, puzzle_lens_1, puzzle_lens_2, puzzle_lens_3, "SAMFEO (5k steps)", "Sampling (softmax, adam)", "log (1 - NED) of Sampling vs. SAMFEO (5k)", "samfeo_5k_log_ned", use_log=True, ned=True)
        # scatter_plot(samfeo_10k_1, samfeo_10k_2, samfeo_10k_3, sampling_soft_1, sampling_soft_2, sampling_soft_3, puzzle_lens_1, puzzle_lens_2, puzzle_lens_3, "SAMFEO (10k steps)", "Sampling (softmax, adam)", "log (1 - NED) of Sampling vs. SAMFEO (10k)", "samfeo_10k_log_ned", use_log=True, ned=True)
        pass

    # bar_plot(sampling_soft, samfeo_5k, samfeo_10k, puzzles, undesignable_ids, "Average NED", "Average NED grouped by Length", "hist_ned.png", ned=True, geom=False, no_undesignable=False)
    # bar_plot(sampling_soft, samfeo_5k, samfeo_10k, puzzles, undesignable_ids, "Average NED", "Average NED grouped by Length (w/o undesignable)", "hist_ned_no_undesignable.png", ned=True, geom=False, no_undesignable=True)
    # bar_plot(sampling_soft, samfeo_5k, samfeo_10k, puzzles, undesignable_ids, "Geometric Average NED", "Geometric Average NED grouped by Length", "hist_ned_geom.png", ned=True, geom=True, no_undesignable=False)
    # bar_plot(sampling_soft, samfeo_5k, samfeo_10k, puzzles, undesignable_ids, "Geometric Average NED", "Geometric Average NED grouped by Length (w/o undesignable)", "hist_ned_geom_no_undesignable.png", ned=True, geom=True, no_undesignable=True)