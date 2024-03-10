import unittest
import subprocess

import sys
import argparse
import numpy as np

import RNA
import threading
import concurrent.futures

from ast import literal_eval
import subprocess

from itertools import product
from collections import defaultdict

LF_PATH = "./LinearFold/linearfold"

def valid_seq(y, targeted):
    """Generate all valid sequences for a given structure y"""
    x = []

    pairs = []
    unpairs = []
    stack = []
    for j, c in enumerate(y):
        if c == '(':
            stack.append(j)
        elif c == ')':
            pairs.append((stack.pop(), j))
        else:
            unpairs.append(j)

    if targeted:
        seq = ['A'] * len(y)
        for combination in product(['CG', 'GC'], repeat=len(pairs)):
            for idx, pair in enumerate(combination):
                i, j = pairs[idx]
                seq[i] = pair[0]
                seq[j] = pair[1]
            x.append("".join(seq))
    else:
        seq = ['A'] * len(y)
        for comb1 in product(['CG', 'GC', 'AU', 'UA', 'GU', 'UG'], repeat=len(pairs)):
            for comb2 in product(['A', 'C', 'G', 'U'], repeat=len(unpairs)):
                for idx, pair in enumerate(comb1):
                    i, j = pairs[idx]
                    seq[i] = pair[0]
                    seq[j] = pair[1]
                
                for idx, c in enumerate(comb2):
                    seq[unpairs[idx]] = c

                x.append("".join(seq))

    return x

def probability(seq, dist):
    nucs_to_idx = {'AA': 0, 'CC': 1, 'GG': 2, 'UU': 3,
                   'CG': 0, 'GC': 1, 'GU': 2, 'UG': 3, 'AU': 4, 'UA': 5}

    p = 1.
    for idx, probs in dist.items():
        p *= probs[nucs_to_idx[seq[idx[0]] + seq[idx[1]]]]

    return p

def expected_energy(rna_struct, init, seed, verbose=False):
    # call main, expected free energy mode
    cmds = f"./main --mode expected_energy --init {init} --seed {seed}"
    rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, input=f"{rna_struct}\n".encode())
    stderr_lines = rt.stderr.decode('utf-8').strip().split('\n')
    stdout_lines = rt.stdout.decode('utf-8').strip().split('\n')
    
    dist = {}
    num_var = rna_struct.count('(') + rna_struct.count('.')
    for line in stdout_lines[1:1+num_var]:
        idx, probs = line.split(': ')
        
        probs = [float(prob) for idx, prob in enumerate(probs.split(' ')) if idx % 2 == 1]

        if idx[0] == '(':
            dist[literal_eval(idx)] = probs[:]
        else:
            dist[(int(idx), int(idx))] = probs[:]
    
    value = 0.
    for line in stderr_lines:
        if line.startswith("Total Energy: "):
            value = float(line.split(": ")[1]) / 100.0

    return value, dist

def expected_energy_brute_force(rna_struct, file, dist, verbose=False):
    with open(f"./deltaG/{file}_eval.txt", 'r') as f:
        lines = f.read().split('\n')
        seqs = np.array([line.split(' ')[0] for line in lines if len(line) > 0])
        seqs_prob = [probability(seq, dist) for seq in seqs]

        deltaG = np.array([float(line.split(' ')[1]) for line in lines if len(line) > 0])
        return np.sum(seqs_prob * deltaG)


# Creating a test case class that inherits from unittest.TestCase
class TestExpectedFreeEnergy(unittest.TestCase):

    # Each test method must start with "test_"
    def test_hp(self):
        rna_struct = "(...)"
        name = "hp"
        init = "uniform"
        seed = "1"
        e1, dist = expected_energy(rna_struct, init, seed)
        e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
        self.assertAlmostEqual(e1, e2, places=3)

    def test_hp1(self):
        rna_struct = ".(...)"
        name = "hp1"
        init = "random"
        seed = "1"
        e1, dist = expected_energy(rna_struct, init, seed)
        e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
        self.assertAlmostEqual(e1, e2, places=3)

    def test_hp2(self):
        rna_struct = "(....)"
        name = "hp2"
        init = "random"
        seed = "1"
        e1, dist = expected_energy(rna_struct, init, seed)
        e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
        self.assertAlmostEqual(e1, e2, places=3)

    def test_hp3(self):
        rna_struct = "(......)"
        name = "hp3"
        init = "random"
        seed = "1"
        e1, dist = expected_energy(rna_struct, init, seed)
        e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
        self.assertAlmostEqual(e1, e2, places=3)

    def test_multi(self):
        rna_struct = "((...)(...))"
        name = "multi"
        init = "uniform"
        seed = "42"
        e1, dist = expected_energy(rna_struct, init, seed)
        e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
        self.assertAlmostEqual(e1, e2, places=4)

    def test_multi_1(self):
        rna_struct = "((...)(...))"
        name = "multi"
        init = "random"
        seed = "1"
        e1, dist = expected_energy(rna_struct, init, seed)
        e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
        self.assertAlmostEqual(e1, e2, places=4)

    def test_multi_2(self):
        rna_struct = "((...)(...))"
        name = "multi"
        init = "random"
        seed = "2"
        e1, dist = expected_energy(rna_struct, init, seed)
        e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
        self.assertAlmostEqual(e1, e2, places=4)
    
    # def test_p8(self):
    #     rna_struct = "((((...))))."
    #     name = "p8"
    #     init = "uniform"
    #     seed = "42"
    #     e1, dist = expected_energy(rna_struct, init, seed)
    #     e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
    #     self.assertAlmostEqual(e1, e2, places=4)

    # def test_p8_1(self):
    #     rna_struct = "((((...))))."
    #     name = "p8"
    #     init = "random"
    #     seed = "1"
    #     e1, dist = expected_energy(rna_struct, init, seed)
    #     e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
    #     self.assertAlmostEqual(e1, e2, places=4)

    # def test_p8_2(self):
    #     rna_struct = "((((...))))."
    #     name = "p8"
    #     init = "random"
    #     seed = "2"
    #     e1, dist = expected_energy(rna_struct, init, seed)
    #     e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
    #     self.assertAlmostEqual(e1, e2, places=4)
    
    # def test_internal(self):
    #     rna_struct = "((.((...)).))"
    #     name = "internal"
    #     init = "uniform"
    #     seed = "42"
    #     e1, dist = expected_energy(rna_struct, init, seed)
    #     e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
    #     self.assertAlmostEqual(e1, e2, places=4)

    # def test_multi2_1(self):
    #     rna_struct = "((...)((...)))"
    #     name = "multi2"
    #     init = "random"
    #     seed = "1"
    #     e1, dist = expected_energy(rna_struct, init, seed)
    #     e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
    #     self.assertAlmostEqual(e1, e2, places=4)

    # def test_multi2_2(self):
    #     rna_struct = "((...)((...)))"
    #     name = "multi2"
    #     init = "random"
    #     seed = "2"
    #     e1, dist = expected_energy(rna_struct, init, seed)
    #     e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
    #     self.assertAlmostEqual(e1, e2, places=4)

    # def test_external_1(self):
    #     rna_struct = "((...))((...))"
    #     name = "external"
    #     init = "random"
    #     seed = "1"
    #     e1, dist = expected_energy(rna_struct, init, seed)
    #     e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
    #     self.assertAlmostEqual(e1, e2, places=4)

    # def test_external_2(self):
    #     rna_struct = "((...))((...))"
    #     name = "external"
    #     init = "random"
    #     seed = "2"
    #     e1, dist = expected_energy(rna_struct, init, seed)
    #     e2 = expected_energy_brute_force(rna_struct, f"{name}", dist)
    #     self.assertAlmostEqual(e1, e2, places=4)

# Running the tests
if __name__ == '__main__':
    unittest.main(argv=[''], defaultTest='TestExpectedFreeEnergy')