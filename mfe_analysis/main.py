import numpy as np
import RNA

NEMO_PREFIX = '../nemo/eterna100_nemo_rp'

puzzle_ids = []
structs = []
with open('../data/eterna/eterna100.txt', 'r') as data:
    lines = data.read().split('\n')
    puzzle_ids = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
    structs = [line.split(' ')[1] for line in lines if len(line) > 0]

nemo_seqs = [[] for _ in range(100)]

with open(f'{NEMO_PREFIX}_{0}.csv', 'r') as data:
    lines = data.read().split('\n')
    nemo_structs = [line.split(',')[1] for line in lines[1:] if len(line) > 0]

for i in range(5):
    with open(f'{NEMO_PREFIX}_{i}.csv', 'r') as data:
        lines = data.read().split('\n')
        seqs = [line.split(',')[2] for line in lines[1:] if len(line) > 0]

        for idx, seq in enumerate(seqs):
            nemo_seqs[idx].append(seq)

# reorder nemo
temp = nemo_seqs[:]
nemo_seqs[95] = temp[98]
nemo_seqs[96] = temp[95]
nemo_seqs[97] = temp[99]
nemo_seqs[98] = temp[97]
nemo_seqs[99] = temp[96]

# for i in range(100):
#     if nemo_structs[i] != structs[i]:
#         print(i, nemo_structs.index(structs[i]))

def print_subopt_result(structure, energy, data):
    ss_list = []
    if not structure == None:
        data['ss_list'].append((energy, structure))
        data['counter'] = data['counter'] + 1

def eval_seq(seq, ss, scale=True):
    subopt_data = { 'counter' : 0, 'sequence' : seq, 'ss_list': []}

    fc = RNA.fold_compound(seq)
    fc.subopt_cb(0, print_subopt_result, subopt_data)

    mfe_structs = [st for e, st in subopt_data['ss_list']]
    is_mfe = ss in mfe_structs
    is_umfe = is_mfe and subopt_data['counter'] == 1

    return seq, is_mfe, is_umfe

nemo_mfe = [False for _ in range(100)]
sampling_mfe = [False for _ in range(100)]
samfeo_mfe = [False for _ in range(100)]

nemo_umfe = [False for _ in range(100)]
sampling_umfe = [False for _ in range(100)]
samfeo_umfe = [False for _ in range(100)]

# seq = 'AAAAAAAAUAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAAAAAGAAAAAAGAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAGAAAAAAGAAAAAAAAAAAAAAAAAAGAGAGACAAAAAAAAAAAAAAAAAAGAAAAAAAAAGAGAAAAAAAAGAAAAAAAAAAAAAAAAAAAGAAAGAAAUAAAAAAGAAGAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAUAAAAAAAAAAAAAAAAAA'
# st  = '................................................................................................................................................................................................................................................................................................................................................................................................................'
# print(eval_seq(seq, st))

def prob(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    pr = fc.pr_structure(ss)
    return pr

def ensemble_defect(seq, ss, scale=True):
    fc = RNA.fold_compound(seq)
    if scale:
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
    fc.pf()
    fc.bpp()
    ed = fc.ensemble_defect(ss)
    return ed

nemo_pyx_seq = ["" for i in range(100)]
nemo_ned_seq = ["" for i in range(100)]

for idx, p_id in enumerate(puzzle_ids):
    puzzle = structs[idx]
    nemo = nemo_seqs[idx]

    best_ned = 10
    best_pyx = -1
    for seq in nemo:
        pyx = prob(seq, puzzle)
        ned = ensemble_defect(seq, puzzle)

        if pyx > best_pyx:
            nemo_pyx_seq[idx] = seq

        if ned < best_ned:
            nemo_ned_seq[idx] = seq


        # _, is_mfe, is_umfe = eval_seq(seq, puzzle)

        # if is_mfe:
        #     nemo_mfe[idx] = True

        # if is_umfe:
        #     nemo_umfe[idx] = True


for x in nemo_pyx_seq:
    print(x)
print()
for x in nemo_ned_seq:
    print(x)

exit(0)

with open('samfeo_mfe.txt', 'r') as f:
    lines = f.read().split('\n')
    samfeo_mfe = [line == 'is_mfe' for line in lines]

with open('sampling_mfe.txt', 'r') as f:
    lines = f.read().split('\n')
    sampling_mfe = [line == 'is_mfe' for line in lines]

with open('samfeo_umfe.txt', 'r') as f:
    lines = f.read().split('\n')
    samfeo_umfe = [line == 'is_umfe' for line in lines]

with open('sampling_umfe.txt', 'r') as f:
    lines = f.read().split('\n')
    sampling_umfe = [line == 'is_umfe' for line in lines]

print("nemo (mfe): ", nemo_mfe.count(True))
print("sampling (mfe): ", sampling_mfe.count(True))
print("samfeo (mfe): ", samfeo_mfe.count(True), '\n')

print("nemo (umfe): ", nemo_umfe.count(True))
print("sampling (umfe): ", sampling_umfe.count(True))
print("samfeo (umfe): ", samfeo_umfe.count(True), '\n')

def intersection(nemo, samfeo, sampling):
    # nemo, samfeo, sampling, nemo & samfeo, samfeo & sampling, nemo & sampling, all
    total = 0
    section = [0 for _ in range(7)]
    for i in range(100):
        if nemo[i] or samfeo[i] or sampling[i]:
            total += 1

        if nemo[i] and samfeo[i] and sampling[i]:
            section[6] += 1
        elif nemo[i] and sampling[i]:
            print("Puzzle id solved by nemo-sampling: ", puzzle_ids[i])
            section[5] += 1
        elif samfeo[i] and sampling[i]:
            print("Puzzle id solved by samfeo-sampling: ", puzzle_ids[i])
            section[4] += 1
        elif nemo[i] and samfeo[i]:
            print("Puzzle id solved by nemo-samfeo: ", puzzle_ids[i])
            section[3] += 1
        elif sampling[i]:
            print("Puzzle id solved by sampling only: ", puzzle_ids[i])
            section[2] += 1
        elif samfeo[i]:
            print("Puzzle id solved by samfeo only: ", puzzle_ids[i])
            section[1] += 1
        elif nemo[i]:
            print("Puzzle id solved by nemo only: ", puzzle_ids[i])
            section[0] += 1

    print("Union: ", total)
    print("All: ", section[6])
    print("Nemo & Sampling: ", section[5])
    print("Samfeo & Sampling: ", section[4])
    print("Nemo & Samfeo: ", section[3])
    print("Sampling: ", section[2])
    print("Samfeo: ", section[1])
    print("Nemo: ", section[0])
    print("")

print("mfe intersections")
intersection(nemo_mfe, samfeo_mfe, sampling_mfe)

print("umfe intersections")
intersection(nemo_umfe, samfeo_umfe, sampling_umfe)
