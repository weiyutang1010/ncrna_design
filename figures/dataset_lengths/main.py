import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb}')

files = ['eterna100_v1.txt', 'rfam27.txt', 'rnasolo.txt']
names = ['eterna100_v1', 'rfam27', 'rnasolo']
datasets = {}

# Get structure lengths from each dataset
for file, name in zip(files, names):
    with open(f'./{file}', 'r') as f:
        lines = f.read().split('\n')
        lines = [line.split(' ')[1] for line in lines if len(line) > 0]
        lengths = [len(line) for line in lines]
        datasets[name] = lengths

def compute_shared_bins(dsets, bin_width=50, min_edge=None, max_edge=None):
    """Compute common bin edges for all datasets."""
    all_vals = np.concatenate(list(dsets.values()))
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
    colors = {'eterna100_v1': '#E35235', 'rfam27': '#fcb777', 'rnasolo': '#4573B4'}

    for name, arr in datasets.items():
        bins = compute_shared_bins(dsets, bin_width=bin_width)
        fig, ax = plt.subplots(figsize=(13, 5.50))
        ax.grid(True, alpha=0.5)
        ax.set_axisbelow(True)
        ax.hist(arr, bins=bins, alpha=0.85, color=colors[name], edgecolor="black")
        
        # ax.set_title(summary(name, arr))
        print(summary(name, arr))

        ax.tick_params(axis='both', which='both', labelbottom=True, labelsize=24)
        
        # only label the x-axis for the last plot
        if name == 'rnasolo':
            ax.set_xlabel("Puzzle Length ($nt$)", fontsize=26)
        ax.set_ylabel("Count", fontsize=26)

        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim(bins[0], bins[-1])

        if name == 'eterna100_v1':
            ax.set_yticks([2, 4, 6, 8, 10, 12, 14, 16, 18])
        if name == 'rfam27':
            ax.set_yticks([1, 2, 3])
        if name == 'rnasolo':
            ax.set_yticks([20, 40, 60, 80, 100, 120, 140])

        fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.97)
        plt.savefig(f'./plots/{name}_lengths.pdf')
        print(f"Figure saved to ./plots/{name}_lengths.pdf")

if __name__ == '__main__':
    plot_faceted_histograms(datasets, bin_width=10, xlim=(0, 400))