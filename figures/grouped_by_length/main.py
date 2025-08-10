#!/usr/bin/env python3

import os, sys
import numpy as np
from scipy.stats import gmean
import matplotlib.pyplot as plt
from matplotlib import collections as matcoll

# for boldsymbol
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

def plot(samplingdesign, samfeo, obj="arith"):
    lines_1 = []
    lines_2 = []
    annotations = []
    for i in range(0, 100, 10):
        start, end = i, i + 10
        # arith mean p(y* | x), geometric mean p(y* | x) excluding undesignables, ned, structural dist, ddG
        obj_idx = {"arith": 2, "geom": 2, "ned": 3, "dist": 4, "ddg": 5}
        idx = obj_idx[obj]

        # geom. mean
        lengths = [samplingdesign[x][0] for x in range(start, end)]
        if obj == "geom":
            samplingdesign_subset = [samplingdesign[x][idx] for x in range(start, end) if samplingdesign[x][1] not in undesignable_ids]
            samfeo_subset = [samfeo[x][idx] for x in range(start, end) if samfeo[x][1] not in undesignable_ids]
            mean_samplingdesign = gmean(samplingdesign_subset)
            mean_samfeo = gmean(samfeo_subset)
        else:
            samplingdesign_subset = [samplingdesign[x][idx] for x in range(start, end)]
            samfeo_subset = [samfeo[x][idx] for x in range(start, end)]
            mean_samplingdesign = np.mean(samplingdesign_subset)
            mean_samfeo = np.mean(samfeo_subset)

        pair=[(lengths[0], mean_samplingdesign), (lengths[-1], mean_samplingdesign)]
        lines_1.append(pair)
        pair=[(lengths[0], mean_samfeo), (lengths[-1], mean_samfeo)]
        lines_2.append(pair)

        dx, dy = 0, 0
        arrow= False

        # arith, ddg, dist
        dx1, dy1, dx2, dy2 = 0, 0, 0, 0
        arrow1, arrow2 = False, False
        if obj == "arith" or obj == "ddg" or obj == "dist":
            if end <= 60:
                continue

            dx1, dy1, dx2, dy2 = 0, 0, 0, 0
            arrow1, arrow2 = False, False
        
        # ned
        if obj == "ned":
            if start == 40 or start == 50:
                continue

            if start == 30:
                arrow1 = True
                dx1, dy1 = -6, 0.014

        # geom. mean
        if obj == "geom":
            if start == 50:
                arrow1, arrow2 = True, True
                dx1, dy1 = -17, 0.07
                dx2, dy2 = -2, 0.1


        if obj == 'ned':
            samplingdesign_mean_label = f'{mean_samplingdesign:.3f}'
            samfeo_mean_label = f'{mean_samfeo:.3f}'
        else:
            samplingdesign_mean_label = f'{mean_samplingdesign:.2f}'
            samfeo_mean_label = f'{mean_samfeo:.2f}'

        valign = 'bottom' if mean_samplingdesign >= mean_samfeo else 'top'
        annotation = ((lengths[0] + lengths[-1]) / 2, mean_samplingdesign, samplingdesign_mean_label, valign, 'blue', arrow1, dx1, dy1)
        annotations.append(annotation)

        valign = 'bottom' if mean_samplingdesign < mean_samfeo else 'top'
        annotation = ((lengths[0] + lengths[-1]) / 2, mean_samfeo, samfeo_mean_label, valign, 'red', arrow2, dx2, dy2)
        annotations.append(annotation)

    if obj == "ned" or obj == "geom":
        plt.rcParams["figure.figsize"] = [10, 6.4] # ned, geom
    else:
        plt.rcParams["figure.figsize"] = [11, 5.50] # others

    linecoll_1 = matcoll.LineCollection(lines_1, color='blue', label='SamplingDesign')
    linecoll_2 = matcoll.LineCollection(lines_2, color='red', label='SAMFEO')
    
    fig, ax = plt.subplots()
    ax.add_collection(linecoll_1)
    ax.add_collection(linecoll_2)

    # annotation
    for x, y, text, valign, color, arrow, dx, dy in annotations:
        if arrow:
            dy *= -1 if valign == 'top' else 1
            ax.annotate(text,
                        xy=(x, y), xycoords='data',
                        xytext=(x + dx, y + dy), textcoords='data',
                        color=color, fontsize=20,
                        arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=color, lw=1.5))
        else:
            if valign == 'top':
                if obj == "arith" or obj == "geom":
                    y -= 0.01 # p(y* | x)
                elif obj == "ned":
                    y -= 0.001 # ned
                else:
                    y -= 0.1 # ddg, dist
            plt.text(x + dx, y + dy, text, horizontalalignment='center', verticalalignment=valign, fontsize=20, color=color)

    # arrow annotation
    for i in range(6, 10):
        line1 = lines_1[i] # sampling line [(x1, y1), (x2, y2)]
        line2 = lines_2[i] # samfeo line

        dy = 0
        # dy = 0.008 # arith.
        # dy = 0.01 # geom
        # dy = 0.0005 # ned

        dy *= -1 if line1[0][1] < line2[0][1] else 1

        # draw line from new_y2 to new_y1
        new_x = (line1[0][0] + line1[1][0]) / 2
        new_y1 = line1[0][1] + dy
        new_y2 = line2[0][1]

        hw, hl = 1.2, 1.4

        threshold = 0.0
        if obj == "arith":
            threshold = 10 # no arrows
        if obj == "geom":
            threshold = 0.02
        if obj == "ned":
            threshold = 0.005
        if obj == "ddg":
            threshold = 1
        if obj == "dist":
            threshold = 0.1

        if abs(new_y2 - new_y1) >= threshold:
            ax.annotate("",
                        xy=(new_x, new_y1), xycoords='data',
                        xytext=(new_x, new_y2), textcoords='data',
                        # color=color, fontsize=17,
                        # arrowprops=dict(arrowstyle=f"->, head_width={hw}, head_length={hl}", connectionstyle="arc3", color="black", lw=5))
                        arrowprops=dict(connectionstyle="arc3", edgecolor="black", facecolor="None", lw=1.5, width=10, headwidth=3*9))

    ax.set_xlabel('Puzzle Length', fontsize=20)

    if obj == "arith" or obj == "geom":
        ax.set_ylabel('$p(\\boldsymbol{y}^\\star \\mid \\boldsymbol{x})$', fontsize=20) # arith and geom
        ax.set_ylim([0.0, 1.0]) # arith
        if obj == "geom":
            ax.set_ylim([-0.08, 1.0]) # geom

    if obj == "ned":
        ax.set_ylabel('NED$(\\boldsymbol{x},\\boldsymbol{y}^\\star)$', fontsize=20)
        ax.set_ylim([0.00, 0.09])

    if obj == "dist":
        ax.set_ylabel('$d(\\text{MFE}(\\boldsymbol{x}),\\boldsymbol{y}^\\star)$', fontsize=20)
        ax.set_ylim([-1.5, 25.0])

    if obj == "ddg":
        ax.set_ylabel('$\\Delta \\Delta G^{\\circ}(\\boldsymbol{x},\\boldsymbol{y}^\\star)$ (kcal/mol)', fontsize=20)
        ax.set_ylim([-0.9, 14.0])

    ax.set_xlim([-5, 415])
    ax.tick_params(labelsize=20)

    if obj in ["ned", "dist", "ddg"]:
        plt.legend(loc='upper left', fontsize=20) # ned, ddg, dist
    else:
        plt.legend(loc='upper right', fontsize=20)
    
    plt.margins(x=0.01, tight=True)
    plt.savefig(f"./plots/{obj}.pdf", bbox_inches="tight") # change
    print(f"Figure saved to ./plots/{obj}.pdf")

# plot(samplingdesign, samfeo, "arith")
# plot(samplingdesign, samfeo, "geom")
plot(samplingdesign, samfeo, "ned")
# plot(samplingdesign, samfeo, "dist")
# plot(samplingdesign, samfeo, "ddg")