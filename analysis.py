import os, sys
import numpy as np
import subprocess

kT = 61.63207755

obj = 1 # 0: - log p(y|x) 1: delta G
if len(sys.argv) > 1:
    obj = int(sys.argv[1])

curr_path = os.path.dirname(os.path.abspath(__file__))
data_path = curr_path + '/data/eterna'
results_path = curr_path + '/results'

def get_prob_full_model(seq, struct):
    cmds = f"./linearfold -V --eval --dangle 2"
    rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, input="\n".join([seq, struct]).encode())
    lines = rt.stdout.decode('utf-8').strip().split('\n')
    delta_G = eval(lines[1].split()[1])

    cmds = f"./linearpartition -V -b 0"
    rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, input=seq.encode())
    lines = rt.stdout.decode('utf-8').strip().split('\n')
    Q = float(lines[1])

    prob = np.exp(delta_G * 100 / -kT) / np.exp(Q * 100 / -kT)
    return prob

with open(data_path + '/short_eterna.txt') as data_file:
    lines = data_file.read().split('\n')
    for line in lines:
        rna_id, rna_struct, sol1, sol2 = line.split(' ')

        if obj == 0:
            result_filepath = results_path + f'/result_{rna_id}.txt'
        elif obj == 1:
            result_filepath = results_path + f'/deltaG_{rna_id}.txt'

        # if file exists

        with open(result_filepath) as result_file:
            result_lines = result_file.read().split('\n')
            final_idx = result_lines.index('Final Sequence')

            final_seq = result_lines[final_idx+2]
            prob_sim = result_lines[final_idx+6] + " (simplified)"
            prob_full = "p(y|x) =   " + str(round(get_prob_full_model(final_seq, rna_struct), 4)) + " (full)"

            print(rna_id)
            print(rna_struct)
            print(final_seq)
            print(prob_sim)
            print(prob_full)
            print()

