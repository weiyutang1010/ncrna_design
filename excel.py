import os, sys
import numpy as np

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
               
                mfe_seq = ""
                if len(r_lines[4].split(':  ')) > 1:
                    mfe_seq = r_lines[4].split(':  ')[1]

                umfe_seq = ""
                if len(r_lines[5].split(':  ')) > 1:
                    umfe_seq = r_lines[5].split(':  ')[1]

                results[rna_id].append((best_pyx, best_pyx_seq, best_pyx_step, best_ned, best_ned_seq, best_ned_step, mfe_seq, umfe_seq, file_path, total_steps))

    total_steps = []
    steps = []
    avg_pyx = []
    avg_ned = []
    mfe_count = 0
    umfe_count = 0
    print("id, name, length, p(y | x), p(y | x) seq, ned, ned seq, is_mfe, mfe seq, is_umfe, umfe seq, init")
    for line in lines:
        rna_id, rna_struct = line[0], line[1]

        if len(results[rna_id]) == 0:
            continue

        result = max(results[rna_id], key=lambda x: x[0])

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

        print(f"{int(rna_id)}", end=",")
        print(f"", end=",")
        print(f"{len(rna_struct)}", end=",")

        print(f"{pyx_step}", end=",")
        print(f"{total_step}", end=",")
        print(f"{pyx:.3f}", end=",")
        print(f"{pyx_seq}", end=",")
        # print(f"{ned_step}", end=",")
        print(f"{ned:.3f}", end=",")
        print(f"{ned_seq}", end=",")

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

    print(f"avg p(y|x) step: {np.average(steps)}")
    print(f"avg total step: {np.average(total_steps)}")
    print(f"average pyx: {np.mean(avg_pyx):.3f}")
    print(f"average ned: {np.mean(avg_ned):.3f}")
    print(f"mfe count: {mfe_count}")
    print(f"umfe count: {umfe_count}")

if __name__ == '__main__':
    main()