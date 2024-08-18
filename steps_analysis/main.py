#!/usr/bin/env python3

import os, sys
import numpy as np
from vienna import *

import RNA

result_folders = [
    # "sampling_pyx_targeted_sm_softmax_eps_075_adam", # gakona
    # "sampling_pyx_targeted_sm_softmax_eps_075_adam_lr_decay", # ironcreek
    "sampling_pyx_uniform_sm_softmax_adam", # beethoven
]

for folder in result_folders:
    folder_path = f'../results/{folder}'

    for file in os.listdir(folder_path):
        file_path = f'{folder_path}/{file}'

        new_file_path = f'./results/{folder}/{file}'
        # if os.path.exists(new_file_path):
        #     continue

        with open(file_path, 'r') as fr:
            with open(new_file_path, 'w') as fw:
                lines = fr.read().split('\n')
                struct = lines[0]

                for idx, line in enumerate(lines):
                    if line.startswith('step: '):
                        obj = line.split(', ')[1].split(': ')[1]
                        integral = line.split(', ')[3].split(': ')[1]
                        best_sample = lines[idx+3].split(' ')[0]

                        
                        integral_result = eval_seq(integral, struct)
                        best_sample_result = eval_seq(best_sample, struct)

                        fw.write(f'{obj}\n')
                        fw.write(f'{best_sample_result[0]} {best_sample_result[1]:.10e} {best_sample_result[2]:.10f} {best_sample_result[3]} {best_sample_result[4]} {best_sample_result[5]} {best_sample_result[6]:.2f}\n')
                        fw.write(f'{integral_result[0]} {integral_result[1]:.10e} {integral_result[2]:.10f} {integral_result[3]} {integral_result[4]} {integral_result[5]} {integral_result[6]:.2f}\n')
                        fw.write('\n')

        print(f"{file_path} done", file=sys.stderr)
