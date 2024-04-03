import os, sys
import numpy as np

lines = []
with open(f'./data/eterna/{sys.argv[2]}.txt', 'r') as file:
    lines = file.read().split('\n')
    lines = [(line.split(' ')[0], line.split(' ')[1]) for line in lines if len(line) > 0]

lines = sorted(lines, key=lambda x: len(x[1]))

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

def main(sample):
    results = {line[0]: [] for line in lines}
    results_seq = {line[0]: [] for line in lines}

    for line in lines:
        rna_id, rna_struct = line[0], line[1]

        x = [-1, -2]
        # file_path_list = [f'./analysis/{sys.argv[1]}_uniform_sm/{rna_id}.txt', f'./analysis/{sys.argv[1]}_targeted_sm/{rna_id}.txt']
        file_path_list = [f'./analysis/{sys.argv[1]}_targeted_sm/{rna_id}.txt']
        for x_seed, file_path in zip(x, file_path_list):
            if not os.path.exists(file_path):
                print(f"{file_path} does not exist")
                continue

            with open(file_path, 'r') as file:
                r_lines = file.read().split('\n')

                if len(r_lines) < 2:
                    continue

                if not sample:
                    best_integral_seq = r_lines[0].split(':  ')[1].split(' ')[0]
                    best_integral_prob = float(r_lines[0].split(':  ')[1].split(' ')[1])
                    results_seq[rna_id].append((best_integral_prob, best_integral_seq, x_seed))
                    results[rna_id].append(best_integral_prob)
                else:
                    best_sample_seq = r_lines[1].split(':  ')[1].split(' ')[0]
                    best_sample_prob = float(r_lines[1].split(':  ')[1].split(' ')[1])
                    results_seq[rna_id].append((best_sample_prob, best_sample_seq, x_seed))
                    results[rna_id].append(best_sample_prob)

    for seed in range(int(sys.argv[3])):
        for line in lines:
            rna_id, rna_struct = line[0], line[1]

            file_path = f'./analysis/{sys.argv[1]}_random_sm_{seed}/{rna_id}.txt'
            if not os.path.exists(file_path):
                print(f"{file_path} does not exist")
                continue

            with open(file_path, 'r') as file:
                r_lines = file.read().split('\n')

                if len(r_lines) < 2:
                    continue

                if not sample:
                    best_integral_seq = r_lines[0].split(':  ')[1].split(' ')[0]
                    best_integral_prob = float(r_lines[0].split(':  ')[1].split(' ')[1])
                    results_seq[rna_id].append((best_integral_prob, best_integral_seq, seed))
                    results[rna_id].append(best_integral_prob)
                else:
                    best_sample_seq = r_lines[1].split(':  ')[1].split(' ')[0]
                    best_sample_prob = float(r_lines[1].split(':  ')[1].split(' ')[1])
                    results_seq[rna_id].append((best_sample_prob, best_sample_seq, seed))
                    results[rna_id].append(best_sample_prob)

    best_seen = []
    print("rna_id, best, avg, var, best_seq")
    for line in lines:
        rna_id, rna_struct = line[0], line[1]
        print(f"{int(rna_id)}", end=",")
        print(f"", end=",")
        print(f"{len(line[1])}", end=",")
        print(f"{np.max(results[rna_id]):.3f}", end=",")
        print(f"{np.mean(results[rna_id]):.3f}", end=",")
        print(f"{np.var(results[rna_id]):.3f}", end=",")
        print(f"{max(results_seq[rna_id], key=lambda x: x[0])[1]}", end=",")
        
        seed = max(results_seq[rna_id], key=lambda x: x[0])[2]
        if seed == -1:
            print(f"uniform", end=",")
        elif seed == -2:
            print(f"targeted", end=",")
        else:
            print(f"random", end=",")
        print(f"{seed}", end="")
        print()
        best_seen.append(np.max(results[rna_id]))


    print(f"average best: {np.mean(best_seen):.3f}\n")

if __name__ == '__main__':
    print("Sample")
    main(True)
    print("Integral")
    main(False)