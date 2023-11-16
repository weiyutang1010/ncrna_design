# Python script for parsing results and obtain p(y | x) under full energy model

import os, sys
import numpy as np
import subprocess

kT = 61.63207755

LF_PATH = "./linearfold"
LP_PATH = "./linearpartition"

def main(path):
    curr_path = os.path.dirname(os.path.abspath(__file__))
    data_path = curr_path + '/data/eterna'
    results_path = curr_path + '/results/' + path + '/'

    def get_prob_full_model(seq, struct):
        cmds = f"{LF_PATH} -V --eval --dangle 2" # path to LinearFold
        rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, input="\n".join([seq, struct]).encode())
        lines = rt.stdout.decode('utf-8').strip().split('\n')
        delta_G = eval(lines[1].split()[1])

        cmds = f"{LP_PATH} -V -b 0" # path to LinearPartition
        rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, input=seq.encode())
        lines = rt.stdout.decode('utf-8').strip().split('\n')
        Q = float(lines[1])

        prob = np.exp(delta_G * 100 / -kT) / np.exp(Q * 100 / -kT)
        return prob

    with open(data_path + '/short_eterna.txt') as data_file:
        lines = data_file.read().split('\n')
        for line in lines:
            rna_id, rna_struct, sol1, sol2 = line.split(' ')
            result_filepath = results_path + f'{rna_id}.txt'

            with open(result_filepath) as result_file:
                result_lines = result_file.read().split('\n')

                final_idx = result_lines.index('Final Sequence') + 2
                final_seq = result_lines[final_idx]
                prob_sim = round(float(result_lines[final_idx+4].split(' ')[4]), 3)
                prob_full = round(get_prob_full_model(final_seq, rna_struct), 3)

                start_idx = result_lines.index('Step, Sequence, -log p(y|x), p(y|x)') + 1
                start_seq = result_lines[start_idx].split(',')[1].strip()
                start_prob_sim = round(float(result_lines[start_idx].split(',')[3].strip()), 3)
                start_prob_full = round(get_prob_full_model(start_seq, rna_struct), 3)

                length = len(final_seq)


                print("Id: ", rna_id)
                print("Length: ", len(rna_struct))
                print(rna_struct)
                print(final_seq)
                print("p(y|x) = ", prob_sim, " (simplified)")
                print("p(y|x) = ", prob_full, " (full)")
                print(rna_struct)
                print(start_seq)
                print("p(y|x) = ", start_prob_sim, " (initial simplified)")
                print("p(y|x) = ", start_prob_full, " (initial full)")
                print()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Missing Result Path.")
        exit()

    result_path = sys.argv[1]
    main(result_path) 
