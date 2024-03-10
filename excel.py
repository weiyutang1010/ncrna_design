import os, sys
import numpy as np

lines = []
with open(f'./data/eterna/eterna_n50.txt', 'r') as file:
    lines = file.read().split('\n')
    lines = [(line.split(' ')[0], line.split(' ')[1]) for line in lines if len(line) > 0]

lines = sorted(lines, key=lambda x: len(x[1]))
results = {line[0]: [] for line in lines}
results_seq = {line[0]: [] for line in lines}

# for line in lines:
#     rna_id = line[0]
#     with open(f'./analysis/{sys.argv[1]}/{rna_id}.txt', 'r') as file:
#         r_lines = file.read().split('\n')

#         best_seq = r_lines[-2].split(': ')[1].split(' ')[1]
#         best_seen = float(r_lines[-2].split(': ')[1].split(' ')[2])
#         init_obj = float(r_lines[-3].split(': ')[1].split(' ')[1])
#         final_obj = float(r_lines[-3].split(': ')[1].split(' ')[2])

#         results_seq[rna_id].append((best_seen, best_seq))
#         results[rna_id].append(best_seen)

for seed in range(30):
    for line in lines:
        rna_id, rna_struct = line[0], line[1]

        file_path = f'./analysis/{sys.argv[1]}_{seed}/{rna_id}.txt'
        if not os.path.exists(file_path):
            print(f"{file_path} does not exist")
            continue

        with open(file_path, 'r') as file:
            r_lines = file.read().split('\n')

            if len(r_lines) < 2:
                continue

            best_seq = r_lines[-2].split(': ')[1].split(' ')[1]
            best_seen = float(r_lines[-2].split(': ')[1].split(' ')[2])
            init_obj = float(r_lines[-3].split(': ')[1].split(' ')[1])
            final_obj = float(r_lines[-3].split(': ')[1].split(' ')[2])

            results_seq[rna_id].append((best_seen, best_seq, seed))
            results[rna_id].append(best_seen)

best_seen = []
# print("rna_id, best, avg, var, best_seq")
for line in lines:
    rna_id, rna_struct = line[0], line[1]
    print(f"{int(rna_id)}", end=",")
    print(f"", end=",")
    print(f"{len(line[1])}", end=",")
    print(f"{np.max(results[rna_id]):.3f}", end=",")
    print(f"{np.mean(results[rna_id]):.3f}", end=",")
    print(f"{np.var(results[rna_id]):.3f}", end=",")
    print(f"{max(results_seq[rna_id], key=lambda x: x[0])[1]}", end=",")
    print(f"{max(results_seq[rna_id], key=lambda x: x[0])[2]}", end="")
    print()
    best_seen.append(np.max(results[rna_id]))
    # print(f"{max(results_seq[rna_id], key=lambda x: x[0])[1]}")
    # print(f"{rna_struct}")


print(f"average best: {np.mean(best_seen):.3f}")