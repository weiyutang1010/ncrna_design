#!/usr/bin/env python3

"""
    Given a stop condition, find the new metrics

    metrics: arith. and geom.^ mean p(y* | x), ned, dist, ddg, # mfe, # umfe

    ^: without undesignable puzzles
"""

import os, sys
import numpy as np
from scipy.stats import gmean
from vienna import *
from collections import defaultdict, deque
import concurrent.futures


import matplotlib.pyplot as plt

import RNA

with open('../data/eterna/eterna100.txt', 'r') as data:
    lines = data.read().split('\n')
    ids = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
    structs = [line.split(' ')[1] for line in lines if len(line) > 0]

undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]

class Result:
    def __init__(self):
        self.best_pyx, self.best_pyx_seq, self.pyx_step_found, self.pyx_file = -0.1, "", -1, ""
        self.best_ned, self.best_ned_seq = 1.1, ""
        self.best_dist, self.best_dist_seq = 1000, ""
        self.best_ddg, self.best_ddg_seq = 100000, ""
        self.mfe_found, self.mfe_seq = False, ""
        self.umfe_found, self.umfe_seq, self.umfe_step = False, "", 0

def update_result(pid, step, sol, file_path, results):
    seq, pyx, ned, is_mfe, is_umfe, dist, ddg = sol.split(' ')
    pyx, ned, dist, ddg = float(pyx), float(ned), int(dist), float(ddg)

    if pyx > results[pid].best_pyx:
        results[pid].best_pyx = pyx
        results[pid].best_pyx_seq = seq
        results[pid].pyx_step_found = step
        results[pid].pyx_file = file_path

    if ned < results[pid].best_ned:
        results[pid].best_ned = ned
        results[pid].best_ned_seq = seq

    if dist < results[pid].best_dist:
        results[pid].best_dist = dist
        results[pid].best_dist_seq = seq

    if ddg < results[pid].best_ddg:
        results[pid].best_ddg = ddg
        results[pid].best_ddg_seq = seq

    if is_mfe == "True":
        results[pid].mfe_found = True
        results[pid].mfe_seq = seq

    if is_umfe == "True" and not results[pid].umfe_found:
        results[pid].umfe_found = True
        results[pid].umfe_seq = seq
        results[pid].umfe_step = step

step_results = {i: defaultdict(Result) for i in range(10001)}

def update_result_2(pid, sol, results):
    if sol.best_pyx > results[pid].best_pyx:
        results[pid].best_pyx = sol.best_pyx

    if sol.best_ned < results[pid].best_ned:
        results[pid].best_ned = sol.best_ned

    if sol.best_ddg < results[pid].best_ddg:
        results[pid].best_ddg = sol.best_ddg


def parse(cutoff, result_folders, stop=False, longest=False):
    results = defaultdict(Result) # {id: (metrics*)}
    steps = {}

    for folder in result_folders:
        folder_path = f'./results/{folder}'

        for file in os.listdir(folder_path):
            file_path = f'{folder_path}/{file}'
            name, _ = os.path.splitext(file)
            if len(name) > 3:
                continue
            pid = int(name)

            # k = 15
            # d = deque()
            # moving_avg = 0.0

            with open(file_path, 'r') as f:
                lines = f.read().split('\n')

                last_sample = [-1.0, -1] # prob, step
                last_ma = [100000, -1] # - log prob, step

                results_2 = defaultdict(Result) # {id: (metrics*)}
                step = 0
                for idx in range(0, len(lines), 4):
                    step = idx // 4
                    if len(lines[idx]) == 0:
                        break


                    obj = float(lines[idx])
                    best_sample = lines[idx+1]
                    # integral = lines[idx+2]

                    length = len(best_sample.split(' ')[0])

                    # record results
                    pyx = float(best_sample.split(' ')[1])
                    update_result(pid, step, best_sample, file_path, results)
                    # update_result(pid, step, integral, file_path, results)

                    update_result(pid, step, best_sample, file_path, results_2)
                    update_result_2(pid, results_2[pid], step_results[step])

                    # if stop condition, break
                    if step >= cutoff:
                        steps[(pid, file_path)] = step
                        break

                    # if pyx == 1.0:
                    #     steps[(pid, file_path)] = step
                    #     break

                    # if stop and (not longest or length >= 200):
                    #     # check new sample found
                    #     if pyx > last_sample[0]:
                    #         last_sample[1] = step
                    #         last_sample[0] = pyx

                    #     # check moving average
                    #     if len(d) >= k:
                    #         moving_avg -= d[0]
                    #         d.popleft()
                    #     moving_avg += obj
                    #     d.append(obj)

                    #     if len(d) == k and moving_avg / k < last_ma[0]:
                    #         last_ma = [moving_avg / k, step]

                    #     # if step >= k and step > last_sample[1] + 100:
                    #     if step >= k and step > last_ma[1] + k:
                    #         steps[(pid, file_path)] = step
                    #         break

                if (pid, file_path) not in steps:
                    steps[(pid, file_path)] = (len(lines) - 1) // 4

                while step <= 8001:
                    update_result_2(pid, results_2[pid], step_results[step])
                    step += 1

                    
    probs = [result.best_pyx for result in results.values()]
    neds = [result.best_ned for result in results.values()]
    dists = [result.best_dist for result in results.values()]
    ddgs = [result.best_ddg for result in results.values()]
    mfe = [result.mfe_found for result in results.values()]
    umfe = [result.umfe_found for result in results.values()]

    probs_no_undesignable = [result.best_pyx for pid, result in results.items() if pid not in undesignable_ids]
    print(np.mean(probs), file=sys.stderr)
    print(gmean(probs_no_undesignable), file=sys.stderr)
    print(np.mean(neds), file=sys.stderr)
    print(np.mean(dists), file=sys.stderr)
    print(np.mean(ddgs), file=sys.stderr)
    print(np.sum(mfe), file=sys.stderr)
    print(np.sum(umfe), file=sys.stderr)

    return results, steps


    # cnt = 0
    # for pid in ids:
    #     if "uniform" in results[pid].pyx_file:
    #         cnt += 1
    # print(cnt)

    # for struct, pid in zip(structs, ids):
    #     result = results[pid]

    #     is_mfe = "is_mfe" if result.mfe_found else ""
    #     mfe_seq = result.mfe_seq

    #     is_umfe = "is_umfe" if result.umfe_found else ""
    #     umfe_seq = result.umfe_seq


        # print(f"{result.best_pyx},{result.best_pyx_seq},{result.best_ned},{result.best_ned_seq},{result.best_dist},{result.best_dist_seq},{result.best_ddg},{result.best_ddg_seq},{is_mfe},{mfe_seq},{is_umfe},{umfe_seq}")

        # with open(result.pyx_file, 'r') as f:
        #     lines = f.read().split('\n')
        #     num_steps = min(2000, (len(lines) - 1) // 4)

        # print(pid, struct, num_steps)

    # print(f"{cutoff} done", file=sys.stderr)

    # return np.mean(probs), gmean(probs_no_undesignable), np.mean(neds)

# check
# result_folders = [
#     "sampling_pyx_targeted_sm_softmax_eps_075_adam",
#     "sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay",
#     "sampling_pyx_uniform_sm_softmax_adam",
# ]
# results_1, steps_1 = parse(2000, result_folders, stop=False)
# # print(f"{results_1[62].umfe_step}", file=sys.stderr)
# print("", file=sys.stderr)

# # result_folders = [
# #     "sampling_pyx_targeted_sm_softmax_eps_075_adam",
# #     "sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay",
# # ]
# results_2, steps_2 = parse(2000, result_folders, stop=True, longest=True)
# print("", file=sys.stderr)

# cnt = 0
# for pid in ids:
#     if results_1[pid].best_pyx > results_2[pid].best_pyx:
#         cnt += 1
#         print(f"{pid} {len(results_1[pid].best_pyx_seq)} {results_1[pid].best_pyx:.4f} {results_2[pid].best_pyx:.4f} {results_1[pid].best_pyx - results_2[pid].best_pyx:.4f}")

# print(cnt)

# for pid in ids:
#     file = results_1[pid].pyx_file
#     best_step_1 = results_1[pid].pyx_step_found
#     total_step_1 = steps_1[(pid, file)]

#     file = results_2[pid].pyx_file
#     best_step_2 = results_2[pid].pyx_step_found
#     total_step_2 = steps_2[(pid, file)]

#     print(pid, best_step_1, total_step_1, best_step_2, total_step_2)

#     if total_step_2 < best_step_1:
#         print(pid, file=sys.stderr)

#     if results_1[pid].umfe_found != results_2[pid].umfe_found:
#         print(f"umfe: {pid}", file=sys.stderr)


# print(results_2[79].best_pyx_seq, file=sys.stderr)

result_folders = [
    "sampling_pyx_targeted_sm_softmax_eps_075_adam",
    "sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay",
    "sampling_pyx_uniform_sm_softmax_adam",
]

results_1, steps_1 = parse(8001, result_folders)
for i in range(1, 8001):
    probs = [result.best_pyx for result in step_results[i].values()]
    probs_no_undesignable = [result.best_pyx for pid, result in step_results[i].items() if pid not in undesignable_ids]
    neds = [result.best_ned for result in step_results[i].values()]
    ddgs = [result.best_ddg for result in step_results[i].values()]
    print(f"{i} {np.mean(probs)} {gmean(probs_no_undesignable)} {np.mean(neds)} {np.mean(ddgs)}")


# for z in range(16):
#     xs = [i for i in range((z*500) + 1, ((z+1)*500) + 1)]
#     arith = []
#     geom = []
#     neds = []

#     with concurrent.futures.ThreadPoolExecutor() as executor:
#         temp = [executor.submit(parse, cutoff) for cutoff in xs]
#         concurrent.futures.wait(temp)

#         for future in temp:
#             temp2 = future.result()
#             arith.append(temp2[0])
#             geom.append(temp2[1])
#             neds.append(temp2[2])
#             ddgs.append(temp2[3])

#     for idx in range(len(xs)):
#         print(f"{xs[idx]} {arith[idx]} {geom[idx]} {neds[idx]} {ddgs[idx]}")

# for i in range(1, 8001):
#     x, y, z = parse(i)
#     arith.append(x)
#     geom.append(y)
#     neds.append(z)
#     print(f"{i} done")



# plt.plot(x, arith)
# plt.plot(x, geom)
# plt.savefig("pyx.pdf")
# plt.clf()

# plt.plot(x, neds)
# plt.savefig("ned.pdf")
