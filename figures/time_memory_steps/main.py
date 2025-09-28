#!/usr/bin/env python3

import os, sys
import numpy as np
import matplotlib.pyplot as plt

# for boldsymbol
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')

with open('time_memory.txt', 'r') as f:
    lines = f.read().split('\n')
    lines = lines[1:] # skip header
    
    lengths = [int(line.split(',')[1]) for line in lines if len(line) > 0]
    time_targeted = [float(line.split(',')[2]) for line in lines if len(line) > 0]
    memory_targeted = [float(line.split(',')[3]) for line in lines if len(line) > 0]
    time_uniform = [float(line.split(',')[4]) for line in lines if len(line) > 0]
    memory_uniform = [float(line.split(',')[5]) for line in lines if len(line) > 0]
    time_both = [float(line.split(',')[6]) for line in lines if len(line) > 0]
    memory_both = [float(line.split(',')[7]) for line in lines if len(line) > 0]

with open('step.txt', 'r') as f:
    lines = f.read().split('\n')
    lines = lines[1:] # skip header

    lengths = [int(line.split(',')[1]) for line in lines if len(line) > 0]
    best_steps = [int(line.split(',')[3]) for line in lines if len(line) > 0]
    total_steps = [int(line.split(',')[4]) for line in lines if len(line) > 0]

def plot_time(lengths, time_targeted, time_uniform, time_both):
    fig, ax = plt.subplots(figsize=(13, 5.50))

    colors = {'orange': '#F47F1E', 'blue': '#1B78B2', 'green': '#2DA248'}

    plt.plot(lengths, np.array(time_both), color=colors['blue'], linestyle='', marker='o', markerfacecolor='None', ms=10, alpha=0.8)
    plt.plot(lengths, np.array(time_uniform), color=colors['green'], linestyle='', marker='o', markerfacecolor='None', ms=10, alpha=0.8)
    plt.plot(lengths, np.array(time_targeted), color=colors['orange'], linestyle='', marker='o', markerfacecolor='None', ms=10, alpha=0.8)

    xl = np.linspace(min(lengths), max(lengths), 100)

    # determine best fit line
    # both initialization
    logx_both, logy_both = np.log(lengths), np.log(time_both)
    z_both = np.polyfit(logx_both, logy_both, 1)
    y_fit_both = np.exp(np.polyval(z_both, np.log(xl)))
    plt.plot(xl, y_fit_both, color=colors['blue'], label=f'curve of best fit (combined): $\\sim |\\boldsymbol{{y}}|^{{{z_both[0]:.2f}}}$')

    # uniform initialization
    logx_uniform, logy_uniform = np.log(lengths), np.log(time_uniform)
    z_uniform = np.polyfit(logx_uniform, logy_uniform, 1)
    y_fit_uniform = np.exp(np.polyval(z_uniform, np.log(xl)))
    plt.plot(xl, y_fit_uniform, color=colors['green'], label=f'curve of best fit (uniform): $\\sim |\\boldsymbol{{y}}|^{{{z_uniform[0]:.2f}}}$')

    # targeted initialization
    logx_targeted, logy_targeted = np.log(lengths), np.log(time_targeted)
    z_targeted = np.polyfit(logx_targeted, logy_targeted, 1)
    y_fit_targeted = np.exp(np.polyval(z_targeted, np.log(xl)))
    plt.plot(xl, y_fit_targeted, color=colors['orange'], label=f'curve of best fit ($\\epsilon$-targeted): $\\sim |\\boldsymbol{{y}}|^{{{z_targeted[0]:.2f}}}$')
    
    # labels
    plt.legend(fontsize=26)
    plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 400])
    plt.xlim(0, 410)
    plt.xlabel('Puzzle Length ($nt$)', fontsize=26)
    plt.ylabel('Time (hour)', fontsize=26)
    plt.tick_params(labelsize=24)

    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.97)
    plt.savefig("./plots/time.pdf")
    print("Figure saved to ./plots/time.pdf")

def plot_memory(lengths, memory_both):
    fig, ax = plt.subplots(figsize=(13, 5.50))

    colors = {'orange': '#F47F1E', 'blue': '#1B78B2', 'green': '#2DA248'}

    plt.plot(lengths, np.array(memory_both), color=colors['blue'], linestyle='', marker='o', markerfacecolor='None', ms=10, alpha=0.8)

    xl = np.linspace(min(lengths), max(lengths), 100)

    # determine best fit line
    # both initialization
    # logx_both, logy_both = np.log(lengths), np.log(memory_both)
    z_both = np.polyfit(lengths, memory_both, 1)
    y_fit_both = np.polyval(z_both, xl)
    plt.plot(xl, y_fit_both, color=colors['blue'])


    # labels
    # plt.legend(fontsize=26)
    plt.yticks([0, 2, 4, 6])
    plt.xticks([0, 50, 100, 150, 200, 250, 300, 350, 400])
    plt.xlim(0, 410)
    plt.xlabel('Puzzle Length ($nt$)', fontsize=26)
    plt.ylabel('Memory (GiB, peak RSS)', fontsize=26)
    plt.tick_params(labelsize=24)

    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.97)
    plt.savefig("./plots/memory.pdf")
    print("Figure saved to ./plots/memory.pdf")

def plot_steps(lengths, best_steps, total_steps):
    fig, ax = plt.subplots(figsize=(13, 5.50))

    x = np.arange(len(lengths))
    w_total = 0.9

    extra_steps = np.array(total_steps) - np.array(best_steps)
    extra_steps = extra_steps.tolist()

    ax.bar(x, extra_steps, color="#FBE594", bottom=best_steps, width=w_total, label="total steps")
    ax.bar(x, best_steps, color="#f08a42", width=w_total, label="best solution found")

    ax.set_xlabel("Puzzle Length ($nt$)", fontsize=26)
    ax.set_ylabel("Steps", fontsize=24)

    step = max(1, len(lengths) // 9)               # ~12 labels max
    ax.set_xticks(x[::step])
    xticklabels = [f'${length}$' for length in lengths[::step]]
    ax.set_xticklabels(xticklabels)

    # nice to read
    # ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(0, max(total_steps))
    ax.legend(loc="lower right", fontsize=24, facecolor='white', framealpha=1)

    plt.margins(x=0.01)
    plt.tick_params(axis='y', labelsize=22)
    plt.tick_params(axis='x', labelsize=26)

    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.97)
    plt.savefig("./plots/steps.pdf")
    print("Figure saved to ./plots/steps.pdf")

if __name__ == '__main__':
    plot_time(lengths, time_targeted, time_uniform, time_both)
    # plot_memory(lengths, memory_both)
    # plot_steps(lengths, best_steps, total_steps)