#!/usr/bin/env python3

import os, sys
import numpy as np
from scipy.stats import gmean
import matplotlib.pyplot as plt
from matplotlib import collections as matcoll

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')


undesignable_ids = [50, 52, 57, 60, 61, 67, 72, 78, 80, 81, 86, 87, 88, 90, 91, 92, 96, 99]
unknown_ids = [68, 97, 100]

file_path_1 = f'./samplingdesign.txt'
file_path_2 = f'./samfeo.txt'

with open(file_path_1, 'r') as f:
    lines = f.read().split('\n')
    lines = [line.split(',') for line in lines[1:] if len(line) > 0]
    # length, puzzle id, p(y*|x), ned, dist, ddG
    samplingdesign = [(int(line[1]), int(line[0]), float(line[3]), float(line[5]), int(line[7]), float(line[9])) for line in lines]
    samplingdesign.sort()

with open(file_path_2, 'r') as f:
    lines = f.read().split('\n')
    lines = [line.split(',') for line in lines[1:] if len(line) > 0]
    # length, puzzle id, p(y*|x), ned, dist, ddG
    samfeo = [(int(line[1]), int(line[0]), float(line[3]), float(line[5]), int(line[7]), float(line[9])) for line in lines]
    samfeo.sort()

def plot(samplingdesign, samfeo, obj):
    obj_idx = {"prob": 2, "prob_log": 2, "ned": 3}
    idx = obj_idx[obj]

    ids = [x[1] for x in samplingdesign]
    lengths = [x[0] for x in samplingdesign]
    samplingdesign = [x[idx] for x in samplingdesign]
    samfeo = [x[idx] for x in samfeo]

    designable_pairs = []
    undesignable_pairs = []
    unknown_pairs = []
    pairs = []

    plt.rcParams["figure.figsize"] = [11, 5.50]
    fig, ax = plt.subplots()

    for i in range(100):
        linestyle = '-'
        if ids[i] in undesignable_ids:
            undesignable_pairs.append([(lengths[i], samplingdesign[i]), (lengths[i], samfeo[i])])
            linestyle = '--'
            color=(0.71,0.4,0.11)
        elif ids[i] in unknown_ids:
            unknown_pairs.append([(lengths[i], samplingdesign[i]), (lengths[i], samfeo[i])])
            linestyle = ':'
            color=(0.5,0.5,0.5)
        else:
            designable_pairs.append([(lengths[i], samplingdesign[i]), (lengths[i], samfeo[i])])
            color='black'

        dy = 0
        if obj == "prob":
            dy = 0.015
        if obj == "ned":
            dy = 0.005
        dy *= -1 if samplingdesign[i] < samfeo[i] else 1

        # from x, y1 to x, y2
        x = lengths[i]
        if obj == "prob_log":
            y1 = samplingdesign[i] * 1e1
        else:
            y1 = samplingdesign[i] + dy
        y2 = samfeo[i]
        
        hw, hl = 0.6, 0.9

        # arrow annotation
        if obj == "prob":
            cond = abs(y2 - y1) >= 0.05
        elif obj == "prob_log":
            cond = y1 / y2 >= 1e2
        elif obj == "ned":
            cond = abs(y2 - y1) >= 0.018


        if cond:
            ax.annotate("",
                        xy=(x, y1), xycoords='data',
                        xytext=(x, y2), textcoords='data',
                        # color=color, fontsize=17,
                        arrowprops=dict(arrowstyle=f"->, head_width={hw}, head_length={hl}", linestyle=linestyle, connectionstyle="arc3", color=color, alpha=0.8, lw=1.8))


# linecoll1 = matcoll.LineCollection(designable_pairs, color='black', linestyle='-', lw=2, alpha=0.5, label='Designable')
# linecoll2 = matcoll.LineCollection(undesignable_pairs, color='black', linestyle='--', lw=2, alpha=0.5, label='Undesignable')
# linecoll3 = matcoll.LineCollection(unknown_pairs, color='black', linestyle=':', lw=2, alpha=0.5, label='Unknown')

# ax.add_collection(linecoll1)
# ax.add_collection(linecoll2)
# ax.add_collection(linecoll3)

    ax.scatter(lengths, samplingdesign, marker='o', label='SamplingDesign', alpha=0.85, s=180, linewidth=2, color='blue', facecolors='None')
    ax.scatter(lengths, samfeo, marker='+', label='SAMFEO', alpha=0.85, s=160, linewidth=2, color='red')
    ax.plot([], [], color='black', linestyle='-', lw=2, alpha=0.85, label='Designable')
    ax.plot([], [], color=(0.71,0.4,0.11), linestyle='--', lw=2, alpha=0.85, label='Undesignable')
    ax.plot([], [], color=(0.5,0.5,0.5), linestyle=':', lw=2, alpha=0.85, label='Unknown')

    # annotations are done in latex
    # annotation (geom.)
    # if obj == "prob_log":
    #     for i in [78, 81, 82, 83, 84, 88, 92, 93, 98]:
    #         dx, dy = -3, 10

    #         if ids[i] == 99:
    #             dx, dy = 23, 1.8e1
    #         if ids[i] == 73:
    #             dx, dy = 36, 0.1
    #         if ids[i] == 100:
    #             dx, dy = 20, 1.8e1
    #         if ids[i] == 91:
    #             dx, dy = 15, 1.8e1
    #         if ids[i] == 90:
    #             dx, dy = 15, 1.8e1
    #         if ids[i] == 76:
    #             dx, dy = 32, 1.5e1

    #         color = 'black'
    #         if ids[i] in undesignable_ids:
    #             color = 'red'
    #         elif ids[i] in unknown_ids:
    #             color = 'blue'
    #         else:
    #             color = 'black'

    #         # 83, 84
    #         ax.annotate(f'\#{ids[i]}',
    #             xy=(lengths[i], samfeo[i] / dy), xycoords='data',
    #             xytext=(dx, 0), textcoords='offset points',
    #             # arrowprops=dict(facecolor='black', shrink=0.05),
    #             horizontalalignment='right', verticalalignment='top', fontsize=15, color=color)

    # # annotation arith.
    # if obj == "prob":
    #     for i in [63, 79, 87]:
    #         dx, dy = 0, 0

    #         if ids[i] == 83:
    #             dx, dy = 35, 0
    #         if ids[i] == 38:
    #             dx, dy = -5, 0
    #         if ids[i] == 74:
    #             dx, dy = -5, 0

    #         color = 'black'
    #         if ids[i] in undesignable_ids:
    #             color = 'red'
    #         elif ids[i] in unknown_ids:
    #             color = 'blue'

    #         # 83, 84
    #         ax.annotate(f'\#{ids[i]}',
    #             xy=(lengths[i], samfeo[i]), xycoords='data',
    #             xytext=(dx, dy), textcoords='offset points',
    #             # arrowprops=dict(facecolor='black', shrink=0.05),
    #             horizontalalignment='right', verticalalignment='top', fontsize=15, color=color, alpha=0.8)

    # # annotation (ned)
    # if obj == "ned":
    #     for i in [78, 82, 84, 93]:
    #         dx, dy = -3, 14

    #         if ids[i] == 73:
    #             dx, dy = 14, 18
    #         if ids[i] == 76:
    #             dx, dy = 16, 18
    #         if ids[i] == 90:
    #             dx, dy = 20, 18


    #         color = 'black'
    #         if ids[i] in undesignable_ids:
    #             color = 'red'
    #         elif ids[i] in unknown_ids:
    #             color = 'blue'

    #         ax.annotate(f'\#{ids[i]}',
    #             xy=(lengths[i], samfeo[i]), xycoords='data',
    #             xytext=(dx, dy), textcoords='offset points',
    #             # arrowprops=dict(facecolor='black', shrink=0.05),
    #             horizontalalignment='right', verticalalignment='top', fontsize=15, color=color, alpha=0.8)

    ax.set_xlim([0, 410])
    ax.set_xlabel('Puzzle Length', fontsize=20)

    plt.grid()
    ax.tick_params(labelsize=20)

    if obj == "prob":
        ax.set_ylabel('$p(\\boldsymbol{y} \\mid \\boldsymbol{x})$', fontsize=20)
        # ax.set_ylabel('$p(\\boldsymbol{y}^\\star \\mid \\boldsymbol{x})$', fontsize=20)
        pass
    elif obj == "prob_log":
        ax.set_ylabel('$p(\\boldsymbol{y} \\mid \\boldsymbol{x})$', fontsize=20)
        # ax.set_ylabel('$p(\\boldsymbol{y}^\\star \\mid \\boldsymbol{x})$', fontsize=20)
        plt.legend(loc='lower right', bbox_to_anchor=(0.8, 0.1))
        ax.set_xlim([280, 405])
        plt.yscale('log')
        ax.set_yticks([1, 1e-8, 1e-16, 1e-24, 1e-32, 1e-40, 1e-48, 1e-56])

        legend = plt.legend(fontsize=19)
        legend.get_texts()[3].set_color((0.71,0.4,0.11))
        legend.get_texts()[4].set_color((0.5,0.5,0.5))
    elif obj == "ned":
        ax.set_yticks([0.00, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35])
        ax.set_ylabel('NED$(\\boldsymbol{x}, \\boldsymbol{y})$', fontsize=19)
        # ax.set_ylabel('NED$(\\boldsymbol{x}, \\boldsymbol{y}^\\star)$', fontsize=19)
        legend = plt.legend(fontsize=19)

        legend.get_texts()[3].set_color((0.71,0.4,0.11))
        legend.get_texts()[4].set_color((0.5,0.5,0.5))

    plt.savefig(f'./plots/{obj}_individual.pdf', bbox_inches='tight')
    print(f'Figure saved to ./plots/{obj}_individual.pdf')

plot(samplingdesign, samfeo, "prob")
plot(samplingdesign, samfeo, "prob_log")
plot(samplingdesign, samfeo, "ned")


