#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

with open("./parse2.txt", 'r') as f:
    lines = f.read().split('\n')

    x = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
    arith = [float(line.split(' ')[1]) for line in lines if len(line) > 0]
    geom = [float(line.split(' ')[2]) for line in lines if len(line) > 0]
    neds = [float(line.split(' ')[3]) for line in lines if len(line) > 0]
    ddgs = [float(line.split(' ')[4]) for line in lines if len(line) > 0]

plt.rcParams["figure.figsize"] = [8, 8]
plt.rcParams["figure.autolayout"] = True
fig, (ax1, ax2, ax3) = plt.subplots(3, 1,  gridspec_kw={'height_ratios': [2, 1, 1]}, sharex=True)

fontsize = 13

ax2.set_xlabel("Step",fontsize=fontsize)
ax1.set_ylabel('Boltzmann Probability',fontsize=fontsize)
# ax1.set_ylim(0.3, 0.61)
ax1.plot(x, arith, label='arith. mean')
ax1.plot(x, geom, label='geom. mean (w/o undesignable puzzles)')
ax1.axvline(x=2000, color='black', linestyle='--', alpha=0.3)

ax2.plot(x, neds, color='green')
ax2.axvline(x=2000, color='black', linestyle='--', alpha=0.3)
ax2.set_ylabel("NED",fontsize=fontsize)
ax2.set_ylim([0, 0.17])

ax3.plot(x, ddgs, color='red')
ax3.axvline(x=2000, color='black', linestyle='--', alpha=0.3)
ax3.set_ylabel(r"Free Energy Gap",fontsize=fontsize)
ax3.set_ylim([0, 10])

ax1.legend(fontsize=fontsize-2)
ax1.set_xlim(-200, 8200)

ax1.tick_params(labelbottom=True, axis='both', which='major', labelsize=12)
ax2.tick_params(labelbottom=False, axis='both', which='major', labelsize=12)
    

plt.title('')


plt.savefig("pyx2.pdf", bbox_inches='tight')
plt.clf()

# plt.xlabel("Step")
# plt.ylabel('Normalized Ensemble Defect')
# plt.axvline(x=2000, color='black', linestyle='--', alpha=0.3)
# plt.plot(x, neds)
# plt.savefig("ned.pdf")