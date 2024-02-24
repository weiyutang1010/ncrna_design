# Python script for parsing results

import os, sys
import numpy as np
import matplotlib.pyplot as plt

import RNA

kT = 61.63207755
curr_path = os.path.dirname(os.path.abspath(__file__))
data_path = curr_path + '/data/eterna/'
results_path = curr_path + '/results'

def get_prob_full_model(seq, struct, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    pr = fc.pr_structure(struct)
    return pr

def get_pyx_full(argv):
    """Get p(y | x) under full energy model"""
    while True:
        rna_seq = input()
        if rna_seq == 'q':
            exit()
        rna_struct = input()
        if rna_struct == 'q':
            exit()

        print("x (sequence): ", rna_seq)
        print("y (structure): ", rna_struct)
        print("p(y | x) = ", get_prob_full_model(rna_seq, rna_struct))

def parse_results(argv):
    """Parse result files"""
    if len(sys.argv) < 4:
        print("Missing Result Path.")
        exit()

    result_foldername = argv[2]
    data_filename = argv[3]
    pyx = 1 if result_foldername[:3] == "pyx" else 0

    with open(data_path + f'/{data_filename}') as data_file:
        lines = data_file.read().split('\n')
        lines = [line.split(' ') for line in lines]

        for line in lines:
            rna_id, rna_struct, sol1, sol2 = line
            result_filepath = f'{results_path}/{result_foldername}/{rna_id}.txt'
            length = len(rna_struct)

            with open(result_filepath) as result_file:
                result_lines = result_file.read().split('\n')

                seqs = []
                frac_obj = []
                objs = []
                integral_objs = []

                score = 0.

                for curr_line in result_lines:
                    if curr_line.startswith("step: "):
                        curr_seq = curr_line.split(', ')[3].split(': ')[1]
                        if len(seqs) == 0 or curr_seq != seqs[-1]:
                            seqs.append(curr_line.split(', ')[3].split(': ')[1])
                            integral_objs.append(float(curr_line.split(', ')[4].split(': ')[1]))


                if len(seqs) == 0:
                    continue

                if pyx:
                    objs = [np.exp(-obj) for obj in objs]
                    integral_objs = [np.exp(-integral) for integral in integral_objs]

                init_simp = integral_objs[0]
                init_full = get_prob_full_model(seqs[0], rna_struct)

                best_rna_seq = ""
                best_prob_sim = 0.0
                best_prob_full = 0.0
                for i, seq in enumerate(seqs):
                    pyx_full = get_prob_full_model(seq, rna_struct)
                    if pyx_full <= 1.0 and pyx_full >= best_prob_full:
                        best_rna_seq = seq
                        best_prob_sim = integral_objs[i]
                        best_prob_full = pyx_full

                print("Id: ", rna_id)
                print("Length: ", len(rna_struct))
                if pyx:
                    print(f'p(y|x) = {init_simp:.3f} (initial simplified)')
                print(f"p(y|x) = {init_full:.3f} (initial full)")
                print(rna_struct)
                print(best_rna_seq)
                if pyx:
                    print(f"p(y|x) = {best_prob_sim:.3f} (final simplified)")
                print(f"p(y|x) = {best_prob_full:.3f} (final full)")
                print()

def get_initial_distributions(argv):
    if len(sys.argv) < 3:
        print("Missing Result Path.")
        exit()

    folder = argv[2]

    with open(data_path + argv[3]) as data_file:
        lines = data_file.read().split('\n')
        lines = [line.split(' ') for line in lines]
        # lines.sort(key=lambda x: int(x[0]))

        # Get initialization
        for line in lines:
            rna_id, rna_struct, sol1, sol2 = line
            result_filepath = f'{results_path}/{folder}/{rna_id}.txt'

            with open(result_filepath) as result_file:
                result_lines = result_file.read().split('\n')
                start_dist_idx = result_lines.index('Starting Distribution')
                length = len(rna_struct)
                for i in range(1, length+1):
                    print(' '.join(result_lines[start_dist_idx + i].split(', ')))

def main(mode, argv):
    funcs = [get_pyx_full, parse_results, get_initial_distributions]
    funcs[mode](argv)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Missing Mode. 0 - get p(y | x), 1 - process results, 2 - get initial distributions")
        exit()
    mode = int(sys.argv[1])

    main(mode, sys.argv) 
