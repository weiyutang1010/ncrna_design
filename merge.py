import os, sys
import numpy as np

from vienna import *

def geo_mean_overflow(iterable):
    return np.exp(np.log(iterable).mean())

lines = []
with open(f'./data/{sys.argv[1]}.txt', 'r') as file:
    lines = file.read().split('\n')
    lines = [(line.split(' ')[0], line.split(' ')[1]) for line in lines if len(line) > 0]

def main():
    results = {line[0]: [] for line in lines}
    no_file_found = True

    # for every puzzle id
    for line in lines:
        rna_id, rna_struct = line[0], line[1]
        file_path_list = [f'./analysis/{sys.argv[i]}/{rna_id}.txt' for i in range(2, len(sys.argv))]

        # iterate through all given folders in ./analysis/ and parse the best solutions
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
                best_pyx, best_pyx_seq   = float(r_lines[4].split(': ')[1].split(' ')[1]), r_lines[4].split(': ')[1].split(' ')[2]
                best_ned, best_ned_seq   = float(r_lines[5].split(': ')[1].split(' ')[1]), r_lines[5].split(': ')[1].split(' ')[2]
                best_dist, best_dist_seq =   int(r_lines[6].split(': ')[1].split(' ')[1]), r_lines[6].split(': ')[1].split(' ')[2]
                best_ddg, best_ddg_seq   = float(r_lines[7].split(': ')[1].split(' ')[1]), r_lines[7].split(': ')[1].split(' ')[2]

                mfe_seq = ""
                num_mfe = int(r_lines[10].split(': ')[1].split(' ')[1])
                if num_mfe > 0:
                    mfe_seq = r_lines[10].split(': ')[1].split(' ')[2]

                umfe_seq = ""
                num_umfe = int(r_lines[11].split(': ')[1].split(' ')[1])
                if num_umfe > 0:
                    umfe_seq = r_lines[11].split(': ')[1].split(' ')[2]

                time = float(r_lines[13].split(': ')[1].split(' ')[1])

                results[rna_id].append((best_pyx, best_pyx_seq, best_ned, best_ned_seq, best_dist, best_dist_seq, best_ddg, best_ddg_seq, mfe_seq, umfe_seq, file_path, total_steps, time))

    if no_file_found:
        print("No file found", file=sys.stderr)
        exit(0)

    total_time = 0.0
    total_steps = []
    steps = []
    avg_pyx = []
    avg_ned = []
    mfe_count = 0
    umfe_count = 0
    avg_dist = []
    avg_ddg = []
    # print("id,length,p(y | x),p(y | x) seq,ned,ned seq,dist,dist seq,ddg,ddg seq,is_mfe,mfe seq,is_umfe,umfe seq")
    print("id,length,p(y|x),ned,dist,ddg,is_mfe,is_umfe")
    for line in lines:
        rna_id, rna_struct = line[0], line[1]

        if len(results[rna_id]) == 0:
            continue

        # take the best solution out of all folders
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
            mfe_seq = result[8] if (result[8] != '') else mfe_seq

        for result in results[rna_id]:
            umfe_seq = result[9] if (result[9] != '') else umfe_seq

        print(f"{int(rna_id):2d}", end=",")
        print(f"{len(rna_struct):3d}", end=", ")

        print(f"{pyx:.3f}", end=", ")
        # print(f"{pyx_seq}", end=",")

        print(f"{ned:.3f}", end=",")
        # print(f"{ned_seq}", end=",")

        print(f"{dist:3d}", end=", ")
        # print(f"{dist_seq}", end=",")

        print(f"{ddg:.2f}", end=",")
        # print(f"{ddg_seq}", end=",")

        if mfe_seq != "":
            print("is_mfe", end=",")
            # print(mfe_seq, end=",")
            mfe_count += 1
        else:
            print(f",", end="")

        if umfe_seq != "":
            print("is_umfe", end=",")
            # print(umfe_seq)
            umfe_count += 1
        else:
            print(f",", end="")
        print()
        
        avg_pyx.append(pyx)
        avg_ned.append(ned)
        avg_dist.append(dist)
        avg_ddg.append(ddg)

        for result in results[rna_id]:
            total_time += result[12]

    undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]
    pyx_no_undesignable = [pyx for pyx, line in zip(avg_pyx, lines) if int(line[0]) not in undesignable_ids]

    print("")
    print("Summary")
    print(f"arith. mean Boltzmann prob.  : {np.mean(avg_pyx):.3f}") # arithmetic mean of Boltzmann probability
    print(f"geom.  mean Boltzmann prob. (w/o undesignable): {geo_mean_overflow(pyx_no_undesignable):.5f}") # geometric mean of Boltzmann probability (excluding undesignable puzzles)
    print(f"average norm. ensemble defect: {np.mean(avg_ned):.4f}") # average normalized ensemble defect
    print(f"average struct. distance     : {np.mean(avg_dist):.4f}") # average structural distance
    print(f"average free energy gap      : {np.mean(avg_ddg):.4f}") # average free energy gap (Delta Delta G)
    print(f"mfe  count: {mfe_count}") # number of puzzles with at least one MFE solution found
    print(f"umfe count: {umfe_count}") # number of puzzles with at least one uMFE solution found
    print(f"total time (s): {total_time:.2f}") # total time taken

if __name__ == '__main__':
    main()