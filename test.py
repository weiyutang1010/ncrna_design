import unittest
import subprocess

import sys
import argparse
import numpy as np

import RNA
import threading
import concurrent.futures

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

def expected_energy(rna_struct, init, verbose=False):
    # call main, expected free energy mode
    cmds = f"./main --mode expected_energy --init {init}"
    rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, input=f"{rna_struct}\n".encode())
    lines = rt.stderr.decode('utf-8').strip().split('\n')
    value = float(lines[-1][14:]) / 100.0

    return value

def expected_energy_brute_force(rna_struct, init, verbose=False):
    targeted = (init == 'targeted')

    cache = defaultdict(lambda: 0.)

    # 1. iterate through all valid sequences
    # 2. compute sum of probability * Delta_G(x, y)
    value = 0.
    seqs = valid_seq(rna_struct, targeted)

    # TODO: parallelize
    for seq in seqs:
        cmds = f"{LF_PATH} -V --eval --verbose"
        rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE, input=f"{seq}\n{rna_struct}".encode())
        lines = rt.stdout.decode('utf-8').strip().split('\n')
        value += float(lines[-1][len(seq)+2:-1])
        # print(seq, 1 / len(seqs), lines[-3])

        for line in lines[:-2]:
            loop = line.split(' : ')[0].split(';')[0][:-3]
            energy = float(line.split(' : ')[1])
            cache[loop] += (1 / len(seqs)) * energy

    # for x in cache:
    #     print(x, cache[x])

    # print(value / len(seqs))

    return value / len(seqs)

# Creating a test case class that inherits from unittest.TestCase
class TestExpectedFreeEnergy(unittest.TestCase):

    # Each test method must start with "test_"
    # def test_n5(self):
    #     rna_struct = "(...)"
    #     init = "uniform"
    #     self.assertAlmostEqual(expected_energy(rna_struct, init), expected_energy_brute_force(rna_struct, init))

    def test_n9_targeted(self):
        rna_struct = "(((...)))"
        init = "targeted"
        self.assertAlmostEqual(expected_energy(rna_struct, init), expected_energy_brute_force(rna_struct, init))

    def test_n12_targeted(self):
        rna_struct = "((((...))))."
        init = "targeted"
        self.assertAlmostEqual(expected_energy(rna_struct, init), expected_energy_brute_force(rna_struct, init))

    def test_n12_multi_targeted(self):
        rna_struct = "((...)(...))"
        init = "targeted"
        self.assertAlmostEqual(expected_energy(rna_struct, init), expected_energy_brute_force(rna_struct, init))

    def test_n10_internal_targeted(self):
        rna_struct = "(.(....).)"
        init = "targeted"
        self.assertAlmostEqual(expected_energy(rna_struct, init), expected_energy_brute_force(rna_struct, init))

    def test_n11_internal_targeted(self):
        rna_struct = "(.(....)..)"
        init = "targeted"
        self.assertAlmostEqual(expected_energy(rna_struct, init), expected_energy_brute_force(rna_struct, init))

    def test_n19_targeted(self):
        rna_struct = "(((....))((...)..))"
        init = "targeted"
        self.assertAlmostEqual(expected_energy(rna_struct, init), expected_energy_brute_force(rna_struct, init))

# Running the tests
if __name__ == '__main__':
    unittest.main(argv=[''], defaultTest='TestExpectedFreeEnergy')