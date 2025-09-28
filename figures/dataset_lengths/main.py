import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')

files = ['eterna100_v1.txt', 'eterna100_v2.txt', 'rfam27.txt', 'rnasolo.txt']
names = ['eterna100', 'eterna100', 'rfam27', 'rnasolo']
datasets = {}

# Get structure lengths from each dataset
for file, name in zip(files, names):
    with open(f'./{file}', 'r') as f:
        lines = f.read().split('\n')
        lines = [line.split(' ')[1] for line in lines if len(line) > 0]
        lengths = [len(line) for line in lines]

        if name not in datasets:
            datasets[name] = [lengths]
        else:
            datasets[name].append(lengths)


def compute_shared_bins(dsets, bin_width=50, min_edge=None, max_edge=None):
    """Compute common bin edges for all datasets."""
    values = []
    for dset, arrs in dsets.items():
        for arr in arrs:
            values.append(arr)

    all_vals = np.concatenate(values)
    if min_edge is None:
        min_edge = max(0, int(np.floor(all_vals.min() / bin_width) * bin_width))
    if max_edge is None:
        max_edge = int(np.ceil(all_vals.max() / bin_width) * bin_width)
    # include the rightmost edge
    bins = np.arange(min_edge, max_edge + bin_width, bin_width)
    return bins

def summary(name, arr):
    """Append N and median to legend labels."""
    arr = np.asarray(arr)
    return f"{name} (n={len(arr)}, median={np.median(arr):.0f})"

def plot_faceted_histograms(dsets, bin_width=50, xlim=None):
    # colors = {'Eterna100': '#E35235', 'Eterna100 (v2)': '#fcb777', 'rfam27': '#4573B4', 'rnasolo': '#2DA248'}
    colors = {'Eterna100': '#E35235', 'Eterna100 (v2)': '#4573B4', 'rfam27': '#fcb777', 'rnasolo': '#2DA248'}

    for name, arrs in datasets.items():
        bins = compute_shared_bins(dsets, bin_width=bin_width)
        fig, ax = plt.subplots(figsize=(13, 5.50))
        ax.grid(True, alpha=0.5)
        ax.set_axisbelow(True)

        if name == 'eterna100':
            ax.hist(arrs[0], bins=bins, alpha=0.2, color=colors['Eterna100'], edgecolor="none")
            ax.hist(arrs[1], bins=bins, alpha=0.2, color=colors['Eterna100 (v2)'], edgecolor="none")
            ax.hist(arrs[0], bins=bins, alpha=0.5, fill=False, edgecolor=colors['Eterna100'], linewidth=3)
            ax.hist(arrs[1], bins=bins, alpha=0.5, fill=False, edgecolor=colors['Eterna100 (v2)'], linewidth=3)
            h1, _ = np.histogram(arrs[0], bins=bins)
            h2, _ = np.histogram(arrs[1], bins=bins)
            ax.stairs(h1, bins, label='Eterna100', color=colors['Eterna100'], linewidth=3)
            ax.stairs(h2, bins, label='Eterna100 (v2)', color=colors['Eterna100 (v2)'], linewidth=3, linestyle='--')
            print(summary(name, arrs[0]))
        else:
            ax.hist(arrs[0], bins=bins, alpha=0.4, color=colors[name], label=name) # fill
            ax.hist(arrs[0], bins=bins, fill=False, edgecolor=colors[name], linewidth=3) # edges
            h1, _ = np.histogram(arrs[0], bins=bins)
            ax.stairs(h1, bins, color=colors[name], linewidth=3)
            print(summary(name, arrs[0]))


        ax.tick_params(axis='both', which='both', labelbottom=True, labelsize=24)
        
        # only label the x-axis for the last plot
        # if name == 'rnasolo':
        ax.set_xlabel("Puzzle Length ($nt$)", fontsize=26)
        ax.set_ylabel("Count", fontsize=26)

        ax.set_xlim([0, 400])

        if name == 'eterna100':
            ax.set_yticks([2, 4, 6, 8, 10, 12, 14, 16, 18])
            plt.legend(fontsize=26)
        if name == 'rfam27':
            ax.set_yticks([1, 2, 3])
        if name == 'rnasolo':
            ax.set_yticks([20, 40, 60, 80, 100, 120, 140])

        fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.97)
        plt.savefig(f'./plots/{name}_lengths.pdf')
        print(f"Figure saved to ./plots/{name}_lengths.pdf")

if __name__ == '__main__':
    plot_faceted_histograms(datasets, bin_width=10, xlim=(0, 400))