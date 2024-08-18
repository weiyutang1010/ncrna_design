#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

with open('../data/eterna/eterna100.txt', 'r') as data:
    lines = data.read().split('\n')
    ids = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
    structs = [line.split(' ')[1] for line in lines if len(line) > 0]

# with open("./step_solved.txt", 'r') as f:
#     lines = f.read().split('\n')

#     x = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
#     best_found = [int(line.split(' ')[1]) for line in lines if len(line) > 0]
#     total_steps = [int(line.split(' ')[2]) for line in lines if len(line) > 0]

#     # Create the figure and axes
#     plt.rcParams["figure.figsize"] = [14, 5.50]
#     fig, ax = plt.subplots()

#     # Plot the bars
#     bar2 = ax.bar(x, total_steps, color='#ffbc00', label='total steps')
#     bar1 = ax.bar(x, best_found, color='#ff5300', label='best solution found')
#     ax.axhline(2000, color='grey', linestyle='--')

#     ax.set_xlabel(r'Puzzles (Sorted by Length$\longrightarrow$)', fontsize=15)
#     ax.set_ylabel('Steps', fontsize=15)
#     # ax.set_title('Stacked Bar Chart')
#     ax.set_xticks([])
#     ax.legend(fontsize=14)
#     ax.tick_params(labelsize=14)
#     ax.set_ylim([0, 9650])

#     plt.margins(x=0.01, tight=True)
#     plt.savefig("steps_1.pdf", bbox_inches="tight")

# comparison

with open("./stop_new.txt", 'r') as f:
    lines = f.read().split('\n')
    lines = [line.split(' ') for line in lines if len(line) > 0]
    lines = [(len(structs[i]), int(line[0]), int(line[1]), int(line[2]), int(line[3]), int(line[4])) for i, line in enumerate(lines)]
    lines.sort()

    # x = [int(line.split(' ')[0]) for line in lines if len(line) > 0]
    x = [i+1 for i in range(100)]
    labels = [len(structs[i]) for i in range(100)]
    best_found_1 = [line[2] for line in lines]
    total_steps_1 = [line[3] for line in lines]
    best_found_2 = [line[4] for line in lines]
    total_steps_2 = [line[5] for line in lines]

    # Create the figure and axes
    plt.rcParams["figure.figsize"] = [14, 5.50]
    fig, ax = plt.subplots()

    # Plot the bars
    bar2 = ax.bar(x, total_steps_2, color='#ffbc00', label='total steps')
    bar1 = ax.bar(x, best_found_2, color='#ff5300', label='best solution found')
    # ax.scatter(x, total_steps_2, color="green", label="new total steps")
    # ax.scatter(x, best_found_2, color="blue", label="new best found")
    # ax.axhline(2000, color='grey', linestyle='--')
    ax.axvline(70.5, color='black', linestyle='--')

    ax.set_xlabel(r'Puzzle Length', fontsize=15)
    ax.set_ylabel('Steps', fontsize=15)
    # ax.set_title('Stacked Bar Chart')
    ax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    ax.set_xticklabels([labels[0], labels[9], labels[19], labels[29], labels[39], labels[49], labels[59], labels[69], labels[79], labels[89], labels[99]])
    # ax.legend(fontsize=14)
    ax.tick_params(labelsize=14)
    ax.set_ylim([0, 2000])

    plt.margins(x=0.01, tight=True)
    plt.savefig("stop_3.pdf", bbox_inches="tight")