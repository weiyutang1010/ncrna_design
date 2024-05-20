import os, sys
import numpy as np
import matplotlib.pyplot as plt

UNIFORM_PARALLEL_FOLDER = '../results/sampling_pyx_uniform_sm_softmax_adam_time_parallel'
UNIFORM_NO_PARALLEL_FOLDER = '../results/sampling_pyx_uniform_sm_softmax_adam_time'
TARGETED_PARALLEL_FOLDER = '../results/sampling_pyx_targeted_sm_softmax_adam_time_parallel'
TARGETED_NO_PARALLEL_FOLDER = '../results/sampling_pyx_targeted_sm_softmax_adam_time'
DATA_FILE = '../data/eterna/eterna100.txt'

puzzle_ids = []
structs = []

with open(DATA_FILE, 'r') as data:
    lines = data.read().split('\n')
    puzzle_ids = [line.split(' ')[0] for line in lines if len(line) > 0]
    structs = [line.split(' ')[1] for line in lines if len(line) > 0]

def parse(folder, label):
    struct_len = []
    time_arr = []
    for p_id, struct in zip(puzzle_ids, structs):
        file_path = f"{folder}/{p_id}.txt"

        if not os.path.exists(file_path):
            continue

        with open(file_path, 'r') as f:
            lines = f.read().split('\n')
            
            time = 0.0
            cnt = 0
            for line in lines:
                if line.startswith('step: '):
                    time += float(line.split(', ')[6].split(': ')[1])
                    cnt += 1

            if cnt == 0:
                continue

            time_avg = time / cnt

            struct_len.append(len(struct))
            time_arr.append(time_avg)
    
    plt.scatter(struct_len, time_arr, marker='x', label=label)

    # best fit of line
    x_values = np.linspace(0, max(struct_len), 1000)
    coefficients = np.polyfit(struct_len, time_arr, 3)
    best_fit = np.polyval(coefficients, x_values)
    plt.plot(x_values, best_fit, linestyle='--')

    print(np.polyval(coefficients, 400) * 2000 / 3600)
    print(f"Label = {label}, Max Time = {max(time_arr)}, Average = {np.mean(time_arr)}, 2000 steps = {max(time_arr) * 2000 / 3600} hr")



parse(UNIFORM_PARALLEL_FOLDER, 'Multi-Threading (56 cores)')
parse(UNIFORM_NO_PARALLEL_FOLDER, 'Single Thread')
# parse(TARGETED_PARALLEL_FOLDER, 'Targeted Distribution (multi-threading)')
# parse(TARGETED_NO_PARALLEL_FOLDER, 'Uniform Distribution')


plt.xlabel('Structure Length')
plt.ylabel('Average Time per Step (sec)')
plt.title('Average Time per Step vs. Structure Length (Sample Size = 2500, Uniform Distribution)')
plt.legend()

save_path = './time.png'
plt.savefig(save_path, format="png", bbox_inches="tight")