#!/usr/bin/env python3

import os, sys
import json
import numpy as np
from scipy.stats import gmean
import matplotlib.pyplot as plt

# for boldsymbol
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')

undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]
unknown_ids = [68, 97, 100]

files = [
    'best_pyx.txt',
    'best_ned.txt',
    'best_ddg.txt',
    'entropy.txt'
]

arith = []
geom = []
neds = []
ddgs = []

def parse_file(file, mean="avg", exclude=False):
    with open(f'./data/{file}', 'r') as f:
        lines = f.read().split('\n')
        if exclude:
            ids = [int(line.split('\t')[0]) for line in lines if len(line) > 0 and int(line.split('\t')[0]) not in undesignable_ids]
            values = [json.loads(line.split('\t')[1]) for line in lines if len(line) > 0 and int(line.split('\t')[0]) not in undesignable_ids]
        else:
            ids = [int(line.split('\t')[0]) for line in lines if len(line) > 0]
            values = [json.loads(line.split('\t')[1]) for line in lines if len(line) > 0]
        
        values = np.array(values)
        values = np.transpose(values)

        if mean == "avg":
            return np.mean(values, axis=1).tolist()
        elif mean == "geom":
            gmean_values = []
            for row in values:
                gmean_values.append(gmean(row))
            return gmean_values
        return []



def plot(steps, arith, geom, neds, ddgs, entropy):
    plt.rcParams["figure.figsize"] = [9, 10]
    plt.rcParams["figure.autolayout"] = True
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1,  gridspec_kw={'height_ratios': [2, 1, 1, 1]}, sharex=True)

    fontsize = 17

    ax4.set_xlabel("Step",fontsize=fontsize)
    ax1.set_ylabel('$p(\\boldsymbol{y}^\\star \\mid \\boldsymbol{x})$',fontsize=fontsize)
    # ax1.set_ylabel('Boltzmann probability',fontsize=fontsize)
    # ax1.set_ylim(0.3, 0.61)
    ax1.plot(steps, arith, label='arith. mean', color='blue')
    ax1.plot(steps, geom, label='geom. mean (w/o undesignable puzzles)', color='darkorange')

    ax1.axhline(0.595, ls="--", lw=1.4, color='blue', alpha=0.8, zorder=0)
    ax1.axhline(0.522, ls="--", lw=1.4, color='darkorange', alpha=0.8, zorder=0)
    # ax1.axvline(x=2000, color='black', linestyle='--', alpha=0.3)

    ax2.plot(steps, neds, color='green')
    ax2.set_ylabel("NED$(\\boldsymbol{x},\\boldsymbol{y}^\\star)$",fontsize=fontsize)
    # ax2.set_ylabel("NED",fontsize=fontsize)
    ax2.axhline(0.035, ls="--", lw=1.4, color='green', alpha=0.8, zorder=0)
    ax2.set_ylim([0, 0.17])

    ax3.plot(steps, ddgs, color='red')
    ax3.set_ylabel("$\\Delta \\Delta G^{\\circ}(\\boldsymbol{x},\\boldsymbol{y}^\\star)$\n(kcal/mol)",fontsize=fontsize)
    # ax3.set_ylabel("Free Energy Gap\n(kcal/mol)",fontsize=fontsize)
    ax3.axhline(0.72, ls="--", lw=1.4, color='red', alpha=0.8, zorder=0)
    ax3.set_ylim([0, 10])

    ax4.plot(steps, entropy, color='dimgrey')
    ax4.set_ylabel("Entropy of $p_{\\boldsymbol{y}^\\star}(\\cdot;\\boldsymbol{\\Theta})$",fontsize=fontsize)
    ax4.axhline(42.9, ls="--", lw=1.6, color='dimgrey', alpha=0.8, zorder=0)
    ax4.set_ylim([0, 190])

    ax1.legend(fontsize=fontsize-2)
    ax1.set_xlim(-25, 2000)

    ax1.tick_params(labelbottom=True, axis='both', which='major', labelsize=17)
    ax2.tick_params(labelbottom=True, axis='both', which='major', labelsize=17)
    ax3.tick_params(labelbottom=True, axis='both', which='major', labelsize=17)
    ax4.tick_params(labelbottom=True, axis='both', which='major', labelsize=17)

    # annotations
    ax1.annotate(
        "$0.595$",
        xy=(1, 0.595), xycoords=("axes fraction", "data"),
        xytext=(6, 0), textcoords="offset points",
        va="center", ha="left",
        fontsize=fontsize,
        color='blue',
        bbox=dict(fc="white", ec="none", alpha=0.6)  # readable over the plot
    )

    ax1.annotate(
        "$0.522$",
        xy=(1, 0.522), xycoords=("axes fraction", "data"),
        xytext=(6, 0), textcoords="offset points",
        va="center", ha="left",
        fontsize=fontsize,
        color='darkorange',
        bbox=dict(fc="white", ec="none", alpha=0.6)  # readable over the plot
    )

    ax2.annotate(
        "$0.035$",
        xy=(1, 0.035), xycoords=("axes fraction", "data"),
        xytext=(6, 0), textcoords="offset points",
        va="center", ha="left",
        fontsize=fontsize,
        color='green',
        bbox=dict(fc="white", ec="none", alpha=0.6)  # readable over the plot
    )

    ax3.annotate(
        "$0.72$",
        xy=(1, 0.72), xycoords=("axes fraction", "data"),
        xytext=(6, 0), textcoords="offset points",
        va="center", ha="left",
        fontsize=fontsize,
        color='red',
        bbox=dict(fc="white", ec="none", alpha=0.6)  # readable over the plot
    )

    ax4.annotate(
        "$42.9$",
        xy=(1, 42.9), xycoords=("axes fraction", "data"),
        xytext=(6, 0), textcoords="offset points",
        va="center", ha="left",
        fontsize=fontsize,
        color='dimgrey',
        bbox=dict(fc="white", ec="none", alpha=0.6)  # readable over the plot
    )

    plt.savefig("./plots/averages.pdf", bbox_inches='tight')
    print("Figure saved to ./plots/averages.pdf")
    plt.cla()

if __name__ == '__main__':
    # compute the average for each step
    steps = np.arange(0, 2000)
    arith = parse_file(files[0])
    geom = parse_file(files[0], mean="geom", exclude=True)
    neds = parse_file(files[1])
    ddgs = parse_file(files[2])
    entropy = parse_file(files[3])

    plot(steps, arith, geom, neds, ddgs, entropy)