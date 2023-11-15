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
            length = len(final_seq)

            start_idx = result_lines.index('Starting Distribution')
            start_dist = result_lines[start_idx+1: start_idx + length + 1]
            start_dist = [temp.split(', ') for temp in start_dist]
            for i in range(length):
                start_dist[i] = [float(x) for x in start_dist[i]]

            temp_idx = result_lines.index('Step, Sequence, -log p(y|x), p(y|x)')
            start_prob_simp = float(result_lines[temp_idx+1].split()[-1])

            start_seq = ["A" for i in range(length)]
            stack = []
            for j, c in enumerate(rna_struct):
                if c == '(':
                    stack.append(j)
                elif c == ')':
                    i = stack.pop()

                    prob_CG = start_dist[i][1] * start_dist[j][2]
                    prob_GC = start_dist[i][2] * start_dist[j][1]

                    if prob_CG > prob_GC:
                        start_seq[i] = 'C'
                        start_seq[j] = 'G'
                    else:
                        start_seq[i] = 'G'
                        start_seq[j] = 'C'

            start_seq = "".join(start_seq)
            start_prob_full = round(get_prob_full_model(start_seq, rna_struct),4)
            # print(start_seq)
            # print(prob)

            prob_sim = result_lines[final_idx+6] + " (simplified)"
            prob_full = "p(y|x) =   " + str(round(get_prob_full_model(final_seq, rna_struct), 4)) + " (full)"

            print("Id: ", rna_id)
            print("Length: ", len(rna_struct))
            print(rna_struct)
            print(final_seq)
            print(prob_sim)
            print(prob_full)
            print(rna_struct)
            print(start_seq)
            print("p(y|x) = ", start_prob_simp, "(initial simplified)")
            print("p(y|x) = ", start_prob_full, "(initial full)")
            print()

