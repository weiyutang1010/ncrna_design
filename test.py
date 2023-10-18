import sys
import copy
import time
import subprocess
import numpy as np

kT = 61.63207755

seq = "GGCGUGGCC"
n = len(seq)
cmds = f"RNAsubopt -e 20 -s"
rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, input=seq.encode())
lines = rt.stdout.decode('utf-8').strip().split('\n')

structs = [line[:len(seq)] for line in lines[1:]]

results = []
for struct in structs:
    cmds = f"./linearfold --eval --dangle 2"
    rt = subprocess.run(cmds.split(), stdout=subprocess.PIPE, input="\n".join([seq, struct]).encode())
    lines = rt.stdout.decode('utf-8').strip().split('\n')
    results.append(lines[1].split())

Q = 0.
scores = [eval(result[1]) for result in results]
for score in scores:
    Q += np.exp((score * 100) / -kT)

results = ["{} \033[92m {}\033[00m".format(result[0], result[1]) for result in results]
print(seq)
print("\n".join(results))
print(f"Free Energy of Ensemble: {(np.log(Q) * -kT) / 100.0:.5f} kcal/mol")
print(f"p(y | x) = {np.exp(scores[0] * 100 / -kT) / Q:.5f}")

