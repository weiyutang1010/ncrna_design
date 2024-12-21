import os, sys
import numpy as np
from scipy.stats import gmean
import matplotlib.pyplot as plt

from vienna import *

lines = []
with open(f'./data/eterna/{sys.argv[1]}.txt', 'r') as file:
    lines = file.read().split('\n')
    lines = [(line.split(' ')[0], line.split(' ')[1]) for line in lines if len(line) > 0]

def main():
    results = {line[0]: [] for line in lines}
    no_file_found = True

    for line in lines:
        rna_id, rna_struct = line[0], line[1]
        file_path_list = [f'./analysis/{sys.argv[i]}/{rna_id}.txt' for i in range(2, len(sys.argv))]

        for file_path in file_path_list:
            if not os.path.exists(file_path):
                print(f"{file_path} does not exist", file=sys.stderr)
                continue

            with open(file_path, 'r') as file:
                r_lines = file.read().split('\n')

                if len(r_lines) < 10:
                    print(f"invalid file read: {file_path}", file=sys.stderr)
                    continue

                no_file_found = False

                total_steps = int(r_lines[0].split(': ')[1].split(' ')[1])
                best_pyx_seq, best_pyx, best_pyx_step = r_lines[3].split(': ')[1].split(' ')[1], float(r_lines[3].split(': ')[1].split(' ')[2]), int(r_lines[3].split(': ')[1].split(' ')[3]) 
                best_ned_seq, best_ned, best_ned_step = r_lines[4].split(': ')[1].split(' ')[1], float(r_lines[4].split(': ')[1].split(' ')[2]), int(r_lines[4].split(': ')[1].split(' ')[3]) 
                best_dist_seq, best_dist, best_dist_step = r_lines[5].split(': ')[1].split(' ')[1], float(r_lines[5].split(': ')[1].split(' ')[2]), int(r_lines[5].split(': ')[1].split(' ')[3]) 
                best_ddg_seq, best_ddg, best_ddg_step = r_lines[6].split(': ')[1].split(' ')[1], float(r_lines[6].split(': ')[1].split(' ')[2]), int(r_lines[6].split(': ')[1].split(' ')[3]) 

                mfe_seq = ""
                num_mfe = int(r_lines[8].split(': ')[1].split(' ')[1])
                if num_mfe > 0:
                    mfe_seq = r_lines[8].split(': ')[1].split(' ')[2]

                umfe_seq = ""
                num_umfe = int(r_lines[8].split(': ')[1].split(' ')[1])
                if num_umfe > 0:
                    umfe_seq = r_lines[8].split(': ')[1].split(' ')[2]

                results[rna_id].append((best_pyx, best_pyx_seq, best_ned, best_ned_seq, best_dist, best_dist_seq, best_ddg, best_ddg_seq, mfe_seq, umfe_seq, file_path, total_steps))

    if no_file_found:
        print("No file found", file=sys.stderr)
        exit(0)

    total_steps = []
    steps = []
    avg_pyx = []
    avg_ned = []
    mfe_count = 0
    umfe_count = 0
    avg_dist = []
    avg_ddg = []
    print("id,length,p(y | x),p(y | x) seq,ned,ned seq,dist,dist seq,ddg,ddg seq,is_mfe,mfe seq,is_umfe,umfe seq")
    for line in lines:
        rna_id, rna_struct = line[0], line[1]

        if len(results[rna_id]) == 0:
            continue

        result = max(results[rna_id], key=lambda x: x[0])
        pyx = result[0]
        pyx_seq = result[1]

        result = min(results[rna_id], key=lambda x: x[2])
        ned = result[2]
        ned_seq = result[3]
        
        result = min(results[rna_id], key=lambda x: x[4])
        dist = result[4]
        dist_seq = result[5]
        
        result = min(results[rna_id], key=lambda x: x[6])
        ddg = result[6]
        ddg_seq = result[7]
        
        mfe_seq, umfe_seq = '', ''
        for result in results[rna_id]:
            mfe_seq = result[8] if (result[8] != '') else ''
            umfe_seq = result[9] if (result[9] != '') else ''

        print(f"{int(rna_id)}", end=",")
        print(f"{len(rna_struct)}", end=",")

        print(f"{pyx:.3f}", end=",")
        print(f"{pyx_seq}", end=",")

        print(f"{ned:.3f}", end=",")
        print(f"{ned_seq}", end=",")

        print(f"{dist}", end=",")
        print(f"{dist_seq}", end=",")

        print(f"{ddg:.2f}", end=",")
        print(f"{ddg_seq}", end=",")

        if mfe_seq != "":
            print("is_mfe", end=",")
            print(mfe_seq, end=",")
            mfe_count += 1
        else:
            print(f",", end=",")

        if umfe_seq != "":
            print("is_umfe", end=",")
            print(umfe_seq)
            umfe_count += 1
        else:
            print(f",")
        
        avg_pyx.append(pyx)
        avg_ned.append(ned)
        avg_dist.append(dist)
        avg_ddg.append(ddg)

    undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]
    pyx_no_undesignable = [pyx for pyx, line in zip(avg_pyx, lines) if int(line[0]) not in undesignable_ids]

    print("")
    print(f"arith. mean pyx: {np.mean(avg_pyx):.3f}", file=sys.stderr)
    print(f"geom. mean pyx (w/o undesignable): {gmean(pyx_no_undesignable):.5f}", file=sys.stderr)
    print(f"average ned: {np.mean(avg_ned):.4f}", file=sys.stderr)
    print(f"average dist: {np.mean(avg_dist):.4f}", file=sys.stderr)
    print(f"average ddg: {np.mean(avg_ddg):.4f}", file=sys.stderr)
    print(f"mfe count: {mfe_count}", file=sys.stderr)
    print(f"umfe count: {umfe_count}", file=sys.stderr)

if __name__ == '__main__':
    main()