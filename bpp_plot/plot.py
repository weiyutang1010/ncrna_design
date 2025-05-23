# python plot.py -x seq.txt -y bpp.txt -z struct.txt --bpp --id 70

import os, sys
import evaluation
import subprocess
from vienna import *

import numpy as np 
from matplotlib import colormaps
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import math
import argparse
from collections import defaultdict
from PIL import Image

plot_data = []

def get_pairs(struct):
    # pairs = [(i_1, j_1), (i_2, j_2), ...] where i_k, j_k is a matching pair in struct
    pairs = []
    stack = []

    for j, c in enumerate(struct):
        if c == '(':
            stack.append(j)
        elif c == ')':
            i = stack.pop()
            pairs.append((i, j))

    return sorted(pairs)

def get_mfe(seq, struct):
    ss_mfe = mfe(seq)[0]
    
    pos = []
    mp, stk = {}, []
    for j, c in enumerate(struct):
        if c == '(':
            stk.append(j)
            mp[j] = -1
        elif c == ')':
            i = stk.pop()
            mp[j] = i
            mp[i] = j
        else:
            mp[j] = j


    stk = []
    for j, c in enumerate(ss_mfe):
        if c == '(':
            stk.append(j)
        elif c == ')':
            i = stk.pop()
            
            if mp[j] == i:
                pos.append(j)
                pos.append(i)
        else:
            if mp[j] == j:
                pos.append(j)

    pos = [x for x in range(len(struct)) if x not in pos]
    pos.sort()
    return pos, ss_mfe


def format_bpp(bpp):
    n = len(bpp)
    bpp_arr = []

    for i, row_i in enumerate(bpp):
        for j, pair_prob in enumerate(row_i):
            if i < j and pair_prob >= 0.01:
                bpp_arr.append((i, j, pair_prob))

    return bpp_arr

def calculate_arc_angles(center, point1, point2):
    # Convert points to angles in radians
    theta1_rad = math.atan2(point1[1] - center[1], point1[0] - center[0])
    theta2_rad = math.atan2(point2[1] - center[1], point2[0] - center[0])
    
    # Normalize angles to [0, 2*pi) range
    theta1_norm = theta1_rad % (2 * math.pi)
    theta2_norm = theta2_rad % (2 * math.pi)
    
    # Convert normalized radians to degrees
    theta1_deg = math.degrees(theta1_norm)
    theta2_deg = math.degrees(theta2_norm)
    
    # Calculate the difference in angles
    diff = theta2_norm - theta1_norm
    
    # Determine if we need to adjust angles to get the smaller arc
    if diff < 0:
        diff += 2 * math.pi  # Adjust if theta2 is before theta1
    
    if diff > math.pi:
        start_deg, end_deg = theta2_deg, theta1_deg
    else:
        start_deg, end_deg = theta1_deg, theta2_deg

    # Ensure angles are in [0, 360) range
    start_deg %= 360
    end_deg %= 360

    return start_deg, end_deg


def calculate_circle_center_and_radius(A, B, C):
    # Transform the input points to numpy arrays for easier manipulation
    A, B, C = np.array(A), np.array(B), np.array(C)
    
    # Using the circumcenter formula
    # The circumcenter is the intersection of the perpendicular bisectors of two sides of the triangle
    d = 2 * (A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]))
    
    # Prevent division by zero in case points are collinear
    if d == 0:
        raise ValueError("Points are collinear")
    
    ux = ((A[0] ** 2 + A[1] ** 2) * (B[1] - C[1]) + (B[0] ** 2 + B[1] ** 2) * (C[1] - A[1]) + (C[0] ** 2 + C[1] ** 2) * (A[1] - B[1])) / d
    uy = ((A[0] ** 2 + A[1] ** 2) * (C[0] - B[0]) + (B[0] ** 2 + B[1] ** 2) * (A[0] - C[0]) + (C[0] ** 2 + C[1] ** 2) * (B[0] - A[0])) / d
    
    center = np.array([ux, uy])
    radius = np.linalg.norm(A - center)
    
    return center, radius

def rotate_point(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.
    
    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

def find_equilateral_points(A, B):
    # Calculate the angle in radians to rotate by 60 degrees (π/3) for one direction
    # and by -60 degrees (-π/3) for the opposite direction to get the two possible C points
    angle_60 = math.pi / 3
    angle_300 = -math.pi / 3  # Equivalent to rotating 60 degrees in the opposite direction
    
    # Rotate B around A by 60 degrees to get one possible position for C
    C1 = rotate_point(A, B, angle_60)
    
    # Rotate B around A by -60 degrees (or 300 degrees) to get the other position for C
    C2 = rotate_point(A, B, angle_300)
    
    return C1, C2

def draw_rna_circle(bpp, seq_len, pairs, folder, puzzle_id, opening_size=0.1, num_index_labels=20, text_offset = 1.025, label_offset = 0.05):
    """
    Draw a circular RNA plot with a customizable opening size, and mark n indices on the circumference based on the sequence size.
    
    Parameters:
    - opening_size: A float representing the fraction of the circle to be left open.
    - sequence_size: The total size of the RNA sequence.
    - n: The number of indices to mark on the circumference, equidistant from each other, following the sequence direction.
    """
    # Calculate the opening angle based on the specified opening size
    opening_angle = 2 * np.pi * opening_size
    
    # Create a new figure for drawing
    plt.figure(figsize=(8, 9))

    # Correctly place labels for 5' and 3' ends
    plt.text(np.sin(opening_angle / 2 - label_offset - 0.05), np.cos(opening_angle / 2 - label_offset - 0.1), "3'", ha='left', va='center', fontsize=13)
    plt.text(-np.sin(opening_angle / 2 - label_offset), np.cos(opening_angle / 2 - label_offset - 0.1), "5'", ha='center', va='center', fontsize=13)
    opening_angle += 0.1

    # Define the angle range for the circle excluding the opening, and for marking indices
    theta = np.linspace(opening_angle / 2, (2 * np.pi) - (opening_angle / 2), 100)
    x = np.sin(theta)
    y = np.cos(theta)
    plt.plot(x, y, color='black')

    # Remove axis lines and ticks
    plt.axis('equal')
    plt.axis('off')

    # Function to get the angle (on the circle) for a given index
    get_angle = lambda index: ((2 * np.pi) - (opening_angle / 2)) - ((index) / (seq_len - 1)) * ((2 * np.pi) - opening_angle)

    for index in np.linspace(0, seq_len - 1, num_index_labels, dtype=int):
        angle = get_angle(index)
        ix = np.sin(angle) * text_offset
        iy = np.cos(angle) * text_offset
        plt.text(ix, iy, str(index + 1), ha='center', va='center', fontsize=11, color='black', rotation=-np.degrees(angle))

    # Draw segments for pairs
    for pair in bpp:
        start_index, end_index, pair_prob = pair[0], pair[1], pair[2]

        # function to get the angle for a given index
        start_angle = get_angle(start_index)
        end_angle = get_angle(end_index)

        # plt.text(np.sin(start_angle), np.cos(start_angle), 'x', ha='center', va='center', fontsize=10, color='blue')
        # plt.text(np.sin(end_angle), np.cos(end_angle), 'x', ha='center', va='center', fontsize=10, color='red')

        start_point = (np.sin(start_angle), np.cos(start_angle))
        end_point = (np.sin(end_angle), np.cos(end_angle))
        mid_points = find_equilateral_points(start_point, end_point)

        # choose a midpoint that lies oustide the circle, if both midpoints lie outside, choose one which is farther from the circle
        if np.linalg.norm(mid_points[0]) > 1 and np.linalg.norm(mid_points[1]) > 1:
            mid_point = max(mid_points, key=lambda x: np.linalg.norm(x))
        else:
            mid_point = mid_points[0] if np.linalg.norm(mid_points[0]) > 1 else mid_points[1]


        def calculate_distance(point1, point2):
            """Calculate the Euclidean distance between two points."""
            return np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

        def calculate_adjustment_factor(dist, diameter):
            """Calculate adjustment factor based on the distance and diameter."""
            # if dist >= diameter:
            #     return 0  # Factor for a straight line 

            if dist >= 0.999 * diameter:
                return (dist ** 6) / diameter
            elif dist >= 0.99 * diameter:
                return (dist ** 5) / diameter
            else:
                return (dist ** 4) / diameter

        def adjust_midpoint(mid_point, factor=0.8):
            direction_vector = np.array(mid_point) / np.linalg.norm(mid_point)  # Unit vector in the direction of center point            
            adjusted_center = mid_point + (direction_vector * factor)  # Move the center point outward
            return adjusted_center
        
        # plt.text(mid_point[0], mid_point[1], 'x', ha='center', va='center', fontsize=10, color='green')

        factor = calculate_adjustment_factor(calculate_distance(start_point, end_point), 2)
        
        plot_data.append((calculate_distance(start_point, end_point), factor))
        mid_point = adjust_midpoint(mid_point, factor)


        # Calculate the radius using the distance to the plot's center
        center, radius = calculate_circle_center_and_radius(start_point, end_point, mid_point)
        theta1, theta2 = calculate_arc_angles(center, start_point, end_point)

        # check ziv range hits
        pair_matched = (start_index, end_index) in pairs

        # Draw the arc
        color = 'blue' if pair_matched else 'red'
        alpha = 1 if not pair_prob else pair_prob
        arc = Arc(center, 2*radius, 2*radius, angle=0, theta1=theta1, theta2=theta2, color=color, alpha=alpha)
        plt.gca().add_patch(arc)
    

    # print title    
    plt.text(0, 1.2, f"Puzzle {puzzle_id}", ha='center', va='center', fontsize=18, color='black')

    # save to file
    save_path = f"./plots/{folder}/circular_plots/{puzzle_id}.png"
    plt.savefig(save_path, dpi=400, bbox_inches="tight")
    print(f"Circular plot saved to {save_path}")

    plt.close()

def draw_rna_linear(bpp, seq_len, pairs, folder, puzzle_id, opening_size=0.1, num_index_labels=10, text_offset = 1.025, label_offset = 0.05, is_mfe=False, mfe_pairs=[]):
    # Create a new figure for drawing
    figwidth, figheight = max(12, seq_len // 15), 10
    plt.figure(figsize=(figwidth, figheight))

    # Correctly place labels for 5' and 3' ends
    label_dist_x, label_dist_y = seq_len * 0.025, 0
    plt.text(0 - label_dist_x, label_dist_y, "5'", ha='left', va='center', fontsize=13)
    plt.text(seq_len-1 + label_dist_x, label_dist_y, "3'", ha='center', va='center', fontsize=13)

    # Remove axis lines and ticks
    plt.axis('equal')
    plt.axis('off')

    # Place indices on the axis
    indices_height = seq_len * (-0.015)
    # plt.text(0, indices_height, str(1), ha='center', va='center', fontsize=13, color='black',)
    # plt.text(seq_len-1, indices_height, str(seq_len), ha='center', va='center', fontsize=13, color='black',)
    
    # index_gap = max(1, seq_len // (num_index_labels-2))
    # for index in range(index_gap, seq_len - 1, index_gap):
    #     plt.text(index, indices_height, str(index + 1), ha='center', va='center', fontsize=13, color='black',)

    for index in np.linspace(0, seq_len - 1, 10).round().astype(int):
        plt.text(index, indices_height, str(index + 1), ha='center', va='center', fontsize=13, color='black',)

    # Draw segments for pairs
    bpp_pairs = []
    for pair in bpp:
        start_index, end_index, pair_prob = pair[0], pair[1], pair[2]
        pair_matched = (start_index, end_index) in pairs
        bpp_pairs.append((start_index, end_index))

        if is_mfe and (start_index, end_index) not in mfe_pairs:
            continue

        # Draw the arc

        # color
        round_prob = math.ceil(pair_prob * 10) / 10
        # colors = {-1.0: "#FF0000", -0.9: "#FF1A1A", -0.8 : "#FF3333", -0.7: "#FF4D4D", -0.6: "#FF6666", -0.5: "#FF8080", -0.4: "#FF9999", -0.3: "#FFB3B3", -0.2: "#FFCCCC", -0.1: "#FFE6E6", 0.0: "#FFA000", 0.1: "#CCCCFF", 0.2: "#B3B3FF", 0.3: "#9999FF", 0.4: "#8080FF", 0.5: "#6666FF", 0.6: "#4D4DFF", 0.7: "#3333FF", 0.8: "#1A1AFF", 0.9: "#0000FF"}


        color = 'blue' if pair_matched else 'red'
        # alpha = 1 if not pair_prob else pair_prob

        # equation for ellipses
        # ref: https://en.wikipedia.org/wiki/Ellipse

        center = ((start_index + end_index) / 2, 0)
        width = (end_index - start_index) / 2
        height = (width+3) * 0.55

        x = np.linspace(start_index, end_index, max(400, seq_len * 4))
        y = np.sqrt((1. - ((x - center[0]) ** 2 / (width * width))) * (height * height)) + center[1]

        if not pair_matched:
            y *= -1

        # plt.plot(x, y, color=color, alpha=alpha)
        # if round_prob == 0.0:
        #     if pair_prob < 0.01:
        #         color = "#FFA000"
        #     else:
        #         color = "#E6E6FF"
        alpha = round_prob
        if pair_matched and pair_prob < 0.1:
            alpha = 1.0
            color = "#FFAA00"
        plt.plot(x, y, alpha=alpha, color=color)

        # if pair_matched and pair_prob >= 0. and pair_prob < .1:
        #     plt.plot(x, y, color="darkorange", alpha=.8)

    
    for pair in pairs:
        if pair not in bpp_pairs:
            start_index, end_index = pair

            center = ((start_index + end_index) / 2, 0)
            width = (end_index - start_index) / 2
            height = (width+3) * 0.55

            x = np.linspace(start_index, end_index, max(400, seq_len * 4))
            y = np.sqrt((1. - ((x - center[0]) ** 2 / (width * width))) * (height * height)) + center[1]

            plt.plot(x, y, color="#FFA000")

    # Draw axis
    x = np.linspace(0, seq_len-1, seq_len)
    y = np.zeros(seq_len)
    plt.plot(x, y, color='black')

    # print title
    # title_height = seq_len * .25
    # plt.text((seq_len - 1) // 2, title_height, f"Puzzle {puzzle_id}", ha='center', va='center', fontsize=18, color='black')
    
    # save to file
    mfe = "_mfe" if is_mfe else ""
    save_format = ".pdf"
    # save_format = ".png"
    save_path = f"./plots/{folder}/linear_plots/{puzzle_id}{mfe}{save_format}"
    plt.savefig(save_path, bbox_inches='tight')
    print(f"Linear plot saved to {save_path}", file=sys.stderr)

    plt.close()

def draw_pos_defects(seq, struct, puzzle_id, defect, folder, layout):
    command = f"echo -e \">{puzzle_id}\\n{seq}\\n{struct}\" | RNAplot -t {layout}"
    command += " --pre \""

    for j, x in enumerate(defect):
        # ref: https://stackoverflow.com/a/7947812
        # white -> red scale (linear interpolation)
        color1 = (255. / 255., 255. / 255., 255. / 255)
        color2 = (255. / 255., 0./ 255., 0. / 255)

        color = (color1[0] + x * (color2[0] - color1[0]), 
                 color1[1] + x * (color2[1] - color1[1]), 
                 color1[2] + x * (color2[2] - color1[2]))

        command += f"{j+1} {j+1} 13 {color[0]} {color[1]} {color[2]} omark "
    
    dx, dy = -1.7, -1.7
    for j in range(9, len(struct), 10):
        command += f"{j+1} {dx} {dy} ({j+1}) Label "

    command += f"1 {-dx - 1.0} {dy + 0.2} (1) Label "

    command += f"{len(struct)} {dx} {dy} ({len(struct)}) Label "

    command += "\""

    save_pdf_path = f"./plots/{folder}/pos_defects/{puzzle_id}.pdf"

    subprocess.call(command, shell=True)
    subprocess.call(f"ps2pdf -dEPSCrop {puzzle_id}_ss.ps", shell=True)
    # subprocess.call(f"pdfcrop {puzzle_id}_ss.pdf", shell=True)
    subprocess.call(f"mv {puzzle_id}_ss.pdf {save_pdf_path}", shell=True)
    subprocess.call(f"rm {puzzle_id}_ss.p*", shell=True)

    print(f"Positional defects plot saved to {save_pdf_path}", file=sys.stderr)

def draw_mfe_plot(seq, struct, puzzle_id, diff, folder, layout, is_mfe=False):
    command = f"echo -e \">{puzzle_id}\\n{seq}\\n{struct}\" | RNAplot -t {layout}"
    command += " --pre \""

    if len(diff) > 0:
        diff_group = []
        start, end = diff[0], diff[0]
        for i in range(1, len(diff)):
            if abs(diff[i] - diff[i-1]) <= 1:
                end = diff[i]
            else:
                diff_group.append((start, end))
                start = diff[i]
                end = diff[i]
        diff_group.append((start, end))
        
        # for x in diff:
        #     command += f"{x+1} {x+1} 13 1.0 0.5 0.5 omark "

        for i, j in diff_group:
            command += f"{i+1} {j+1} 13 1.0 0.5 0.5 omark "
    
    dx, dy = -1.7, -1.7
    for j in range(9, len(struct)-1, 10):
        command += f"{j+1} {dx} {dy} ({j+1}) Label "

    command += f"1 {-dx - 2.4} {dy + 0.4} (1) Label "

    command += f"{len(struct)} {dx + 1.9} {dy + 0.4} ({len(struct)}) Label "

    command += "\""

    is_mfe = "_mfe" if is_mfe else ""
    save_pdf_path = f"./plots/{folder}/mfe/{puzzle_id}{is_mfe}.pdf"

    subprocess.call(command, shell=True)
    subprocess.call(f"ps2pdf -dEPSCrop {puzzle_id}_ss.ps", shell=True)
    # subprocess.call(f"pdfcrop {puzzle_id}_ss.pdf", shell=True)
    subprocess.call(f"mv {puzzle_id}_ss.pdf {save_pdf_path}", shell=True)
    subprocess.call(f"rm {puzzle_id}_ss.p*", shell=True)

    print(f"mfe plot saved to {save_pdf_path}", file=sys.stderr)

def create_folders(folder_name, args):
    if not os.path.exists(f"./plots"):
        print("Created directory \"./plots\"", file=sys.stderr)
        os.makedirs(f"./plots")

    if not os.path.exists(f"./plots/{folder_name}"):
        print(f"Created directory \"./plots/{folder_name}\"", file=sys.stderr)
        os.makedirs(f"./plots/{folder_name}")

    if args.linear and not os.path.exists(f"./plots/{folder_name}/linear_plots"):
        print(f"Created directory \"./plots/{folder_name}/linear_plots\"", file=sys.stderr)
        os.makedirs(f"./plots/{folder_name}/linear_plots")

    if args.circular and not os.path.exists(f"./plots/{folder_name}/circular_plots"):
        print(f"Created directory \"./plots/{folder_name}/circular_plots\"", file=sys.stderr)
        os.makedirs(f"./plots/{folder_name}/circular_plots")

    if args.positional and not os.path.exists(f"./plots/{folder_name}/pos_defects"):
        print(f"Created directory \"./plots/{folder_name}/pos_defects\"", file=sys.stderr)
        os.makedirs(f"./plots/{folder_name}/pos_defects")

    if args.mfe and not os.path.exists(f"./plots/{folder_name}/mfe"):
        print(f"Created directory \"./plots/{folder_name}/mfe\"", file=sys.stderr)
        os.makedirs(f"./plots/{folder_name}/mfe")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw RNA plots')
    parser.add_argument('-o', '--folder', default="temp", type=str, help='Output folder')
    parser.add_argument('-c', '--circular', action="store_true", default=False, help='Draw circular RNA plot')
    parser.add_argument('-l', '--linear', action="store_true", default=False, help='Draw linear RNA plot')
    parser.add_argument('-p', '--positional', action="store_true", default=False, help='Draw positional defects RNA plot')
    parser.add_argument('-m', '--mfe', action="store_true", default=False, help='Draw mfe plot')
    parser.add_argument('-t', '--layout', type=str, default='4', help='RNA plot layout')
    
    args = parser.parse_args()

    create_folders(args.folder, args)

    lines = []
    with open('input', 'r') as f:
        lines = f.read().split('\n')
        n = len(lines)

        puzzles_ids = lines[:n//3]
        seqs = lines[n//3:2*n//3]
        structs = lines[2*n//3:]

    for puzzle_id, seq, struct in zip(puzzles_ids, seqs, structs):
        print(f"Puzzle {puzzle_id} start", file=sys.stderr)

        bpp, pos_defect = position_defect(seq, struct)
        # bpp = [(i, j, pair_prob), ...]
        bpp = format_bpp(bpp)
        # pairs = [(i, j), ...]
        pairs = get_pairs(struct)

        diff, mfe_struct = get_mfe(seq, struct)
        mfe_pairs = get_pairs(mfe_struct)

        # print paired prob (w/ incorrect)
        # for idx, x in enumerate(bpp):
        #     i, j, pair_prob = x
        #     if (i, j) in mfe_pairs:
        #         mult = 1
        #         if (i, j) not in pairs:
        #             mult = -1
                
        #     print(i+1, j+1, pair_prob)
        # print()

        # print pos defect
        # for idx, x in enumerate(pos_defect):
        #     print(idx+1, x)
        
        # circular plot
        if args.circular:
            draw_rna_circle(bpp, len(seq), pairs, args.folder, puzzle_id)

        # linear plot
        if args.linear:
            draw_rna_linear(bpp, len(seq), pairs, args.folder, puzzle_id)

            # if args.mfe:
            #     mfe_pairs = get_pairs(mfe_struct)
            #     draw_rna_linear(bpp, len(seq), pairs, args.folder, puzzle_id, is_mfe=True, mfe_pairs=mfe_pairs)

        # positional defects
        if args.positional:
            draw_pos_defects(seq, struct, puzzle_id, pos_defect, args.folder, args.layout, norm=False)

        if args.mfe:
            # sampling_seq = "GGGCAAAGCCCGGCACUGGAAACAGAGGCAAGGGACCCGAAAGGGGAGCGCCCAAAGGGCAACGGGAAACCCGGCUCCGGGAAACCGUCCCAAGAGAGGCGAAAGCGCGCCGCCCUAAGGGCAAGGUCAAAGACCGGCGCUGGAAACACUCUCAACCUCCUGGAAACACUGCCGCCCAAAGGGCUACGGCAAAGCCGGGCAGGCGAAAGCGGAGGAACCGGCGCGAAAGCGGGCCGCGGAAUCCGCAACUCCAAAGGAGGGCCCAGGAAACUGCCGGAAGGAGCCCGAAAGGGAGGCGCCCCAAGGGCAAGGGCAAAGCCCGCCUCGGGAAACCGCUCCAAGCCUCGCGAAAGCGUGCCCCGGAAACCGG"
            # draw_mfe_plot(sampling_seq, struct, puzzle_id, diff, args.folder, args.layout, is_mfe=False)

            draw_mfe_plot(seq, struct, puzzle_id, diff, args.folder, args.layout, is_mfe=False)
            draw_mfe_plot(seq, mfe_struct, puzzle_id, diff, args.folder, args.layout, is_mfe=True)



