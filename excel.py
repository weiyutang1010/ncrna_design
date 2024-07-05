import os, sys
import numpy as np
from scipy.stats import gmean
import matplotlib.pyplot as plt

from vienna import *

lines = []
with open(f'./data/eterna/{sys.argv[1]}.txt', 'r') as file:
    lines = file.read().split('\n')
    lines = [(line.split(' ')[0], line.split(' ')[1]) for line in lines if len(line) > 0]

lines = sorted(lines, key=lambda x: len(x[1]))

def main():
    results = {line[0]: [] for line in lines}

    for line in lines:
        rna_id, rna_struct = line[0], line[1]

        file_path_list = [f'./analysis/{sys.argv[i]}/{rna_id}.txt' for i in range(2, len(sys.argv))]

        for file_path in file_path_list:
            if not os.path.exists(file_path):
                print(f"{file_path} does not exist", file=sys.stderr)
                continue

            with open(file_path, 'r') as file:
                r_lines = file.read().split('\n')

                if len(r_lines) < 6:
                    print(f"invalid file read: {file_path}", file=sys.stderr)
                    continue

                if len(r_lines[2].split(': ')[1].split(' ')) <= 3:
                    continue

                total_steps = int(r_lines[0].split(': ')[1].split(' ')[1])
                best_pyx_seq, best_pyx, best_pyx_step = r_lines[2].split(': ')[1].split(' ')[1], float(r_lines[2].split(': ')[1].split(' ')[2]), int(r_lines[2].split(': ')[1].split(' ')[3]) 
                best_ned_seq, best_ned, best_ned_step = r_lines[3].split(': ')[1].split(' ')[1], float(r_lines[3].split(': ')[1].split(' ')[2]), int(r_lines[3].split(': ')[1].split(' ')[3]) 
                
                # best_dist_seq, best_dist, best_dist_step = r_lines[6].split(': ')[1].split(' ')[1], float(r_lines[6].split(': ')[1].split(' ')[2]), int(r_lines[6].split(': ')[1].split(' ')[3]) 
                # best_ediff_seq, best_ediff, best_ediff_step = r_lines[7].split(': ')[1].split(' ')[1], float(r_lines[7].split(': ')[1].split(' ')[2]), int(r_lines[7].split(': ')[1].split(' ')[3]) 

                best_dist_seq, best_dist, best_dist_step = 0, 0, 0
                best_ediff_seq, best_ediff, best_ediff_step = 0, 0, 0


                mfe_seq = ""
                if len(r_lines[4].split(':  ')) > 1:
                    mfe_seq = r_lines[4].split(':  ')[1]
                    if len(mfe_seq.split(' ')) > 1:
                        mfe_seq = mfe_seq.split(' ')[1]

                    if mfe_seq == '0':
                        mfe_seq = ''

                umfe_seq = ""
                if len(r_lines[5].split(':  ')) > 1:
                    umfe_seq = r_lines[5].split(':  ')[1]
                    if len(umfe_seq.split(' ')) > 1:
                        umfe_seq = umfe_seq.split(' ')[1]

                    if umfe_seq == '0':
                        umfe_seq = ''

                results[rna_id].append((best_pyx, best_pyx_seq, best_pyx_step, best_ned, best_ned_seq, best_ned_step, mfe_seq, umfe_seq, file_path, total_steps, best_dist, best_dist_seq, best_ediff, best_ediff_seq))

    total_steps = []
    steps = []
    avg_pyx = []
    avg_ned = []
    mfe_count = 0
    umfe_count = 0
    avg_dist = []
    avg_ediff = []
    print("id, name, length, p(y | x), p(y | x) seq, ned, ned seq, is_mfe, mfe seq, is_umfe, umfe seq, init")
    for line in lines:
        rna_id, rna_struct = line[0], line[1]

        if len(results[rna_id]) == 0:
            continue

        # result = max(results[rna_id], key=lambda x: x[0]) # for p(y | x)
        result = min(results[rna_id], key=lambda x: x[3]) # for ned

        pyx = result[0]
        pyx_seq = result[1]
        pyx_step = result[2]
        ned = result[3]
        ned_seq = result[4]
        ned_step = result[5]
        mfe_seq = result[6]
        umfe_seq = result[7]
        file_path = result[8]
        total_step = result[9]
        dist = result[10]
        dist_seq = result[11]
        ediff = result[12]
        ediff_seq = result[13]

        print(f"{int(rna_id)}", end=",")
        print(f"", end=",")
        print(f"{len(rna_struct)}", end=",")

        # print(f"{pyx_step}", end=",")
        # print(f"{total_step}", end=",")
        print(f"{pyx:.3f}", end=",")
        print(f"{pyx_seq}", end=",")
        # print(f"{ned_step}", end=",")
        print(f"{ned:.3f}", end=",")
        print(f"{ned_seq}", end=",")

        # print(f"{dist}", end=",")
        # print(f"{dist_seq}", end=",")

        # print(f"{ediff}", end=",")
        # print(f"{ediff_seq}", end=",")

        if mfe_seq != "":
            print("is_mfe", end="")
            mfe_count += 1
        print(f",{mfe_seq}", end=",")

        if umfe_seq != "":
            print("is_umfe", end="")
            umfe_count += 1
        print(f",{umfe_seq}", end=",")

        # check initialization
        if "targeted" in file_path:
            if "050" in file_path:
                print("targeted (eps=0.5)")
            elif "075" in file_path:
                if "lr_decay" in file_path:
                    print("targeted (eps=0.75+lr_decay)")
                else:
                    print("targeted (eps=0.75)")
            else:
                print("targeted")
        elif "uniform" in file_path:
            print("uniform")
        else:
            print("random")
        
        total_steps.append(total_step)
        steps.append(pyx_step)
        avg_pyx.append(pyx)
        avg_ned.append(ned)
        avg_dist.append(dist)
        avg_ediff.append(ediff)

    undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]
    pyx_no_undesignable = [pyx for pyx, line in zip(avg_pyx, lines) if int(line[0]) not in undesignable_ids]

    # print(f"avg p(y|x) step: {np.average(steps)}")
    print(f"avg total step: {np.average(total_steps)}")
    print(f"average pyx: {np.mean(avg_pyx):.3f}")
    print(f"average pyx (geometric w/o undesignable): {gmean(pyx_no_undesignable):.5f}")
    print(f"average ned: {np.mean(avg_ned):.4f}")
    print(f"mfe count: {mfe_count}")
    print(f"umfe count: {umfe_count}")
    print(f"average dist: {np.mean(avg_dist):.4f}")
    print(f"average ediff: {np.mean(avg_ediff):.4f}")

    # plot p(y | x) against NED
    # avg_ned = np.array(avg_ned)
    # plt.plot(np.log(1 - avg_ned), np.log(avg_pyx), marker='x', linestyle='', label='puzzle')
    # plt.xlabel('log(1 - NED)')
    # plt.ylabel('log p(y | x)')
    # plt.title('p(y | x) vs. NED of Eterna100 Solutions')
    # plt.savefig('./pyx_ned.png', dpi=400, bbox_inches='tight')
    # plt.close()

if __name__ == '__main__':
    main()