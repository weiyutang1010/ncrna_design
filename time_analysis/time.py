import os, sys
import numpy as np
import matplotlib.pyplot as plt

UNIFORM_PARALLEL_FOLDER = '../results/sampling_pyx_uniform_sm_softmax_adam_time_parallel'
UNIFORM_NO_PARALLEL_FOLDER = '../results/sampling_pyx_uniform_sm_softmax_adam_time'
UNIFORM_PARALLEL_CACHE_FOLDER = '../results/sampling_pyx_uniform_sm_softmax_adam_time_parallel_cache'
UNIFORM_NED_PARALLEL_FOLDER = '../results/sampling_log_ned_uniform_sm_softmax_adam_time'
UNIFORM_NED_LAZY_PARALLEL_FOLDER = '../results/sampling_log_ned_uniform_sm_softmax_adam_lazy_time'
TARGETED_PARALLEL_FOLDER = '../results/sampling_pyx_targeted_sm_softmax_adam_time_parallel'
TARGETED_NO_PARALLEL_FOLDER = '../results/sampling_pyx_targeted_sm_softmax_adam_time'
TARGETED_PARALLEL_CACHE_FOLDER = '../results/sampling_pyx_targeted_sm_softmax_adam_time_parallel_cache'
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

    struct_len = struct_len
    time_arr = time_arr

    print(len(time_arr))

    plt.scatter(struct_len, time_arr, marker='x', label=label)

    # best fit of line
    x_values = np.linspace(min(struct_len), max(struct_len), 1000)
    coefficients = np.polyfit(struct_len, time_arr, 3)
    best_fit = np.polyval(coefficients, x_values)
    plt.plot(x_values, best_fit, linestyle='--')

    print(np.polyval(coefficients, max(struct_len)))
    print(f"Label = {label}, Max Time = {max(time_arr)}, Average = {np.mean(time_arr)}, 2000 steps = {max(time_arr) * 2000 / 3600} hr")

plt.figure(figsize=(12,6))


parse(UNIFORM_PARALLEL_FOLDER, 'obj = p(y | x)')
parse(UNIFORM_NED_PARALLEL_FOLDER, 'obj = NED')
parse(UNIFORM_NED_LAZY_PARALLEL_FOLDER, 'obj = NED (lazyoutside)')
# parse(UNIFORM_PARALLEL_CACHE_FOLDER, 'Multi-Threading (28 cores) + Caching')
# parse(UNIFORM_NO_PARALLEL_FOLDER, 'Single Thread')
# parse(TARGETED_PARALLEL_CACHE_FOLDER, 'Multi-Threading (28 cores) + Cache')
# parse(TARGETED_PARALLEL_FOLDER, 'Multi-Threading (28 cores)')
# parse(TARGETED_NO_PARALLEL_FOLDER, 'Uniform Distribution')

plt.xlabel('Structure Length')
plt.ylabel('Time per Step (sec)')
plt.title('Avg Time per Step vs. Length (Sample Size=2500, 28 cores)')
plt.legend()

save_path = './time.png'
plt.savefig(save_path, format="png", bbox_inches="tight")
