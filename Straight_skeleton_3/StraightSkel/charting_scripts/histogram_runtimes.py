import os
import matplotlib.pyplot as plt
import argparse
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import matplotlib.ticker as mticker


def plot_runtime_histogram(root_dir):
    runtimes = []
    for subdir_name in os.listdir(root_dir):
        subdir_path = os.path.join(root_dir, subdir_name)

        if not os.path.isdir(subdir_path):
            continue

        runtime_file_path = os.path.join(subdir_path, 'runtime.txt')

        if os.path.exists(runtime_file_path):
            try:
                with open(runtime_file_path, 'r') as runtime_file:
                    runtime = float(runtime_file.readline().strip())
                runtimes.append(runtime)
            except ValueError as e:
                print(f"Error processing file in directory {subdir_name}: {e}")
                continue

    if runtimes:
        max_runtime = max(runtimes)
        num_bins = 20
        bins = np.linspace(0, max_runtime, num_bins + 1)

        hist, _ = np.histogram(runtimes, bins=bins)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        fig, ax = plt.subplots(figsize=(10, 6))  # Adjust figure size

        cmap = LinearSegmentedColormap.from_list("GreenYellowRed", ["green", "yellow", "red"], N=256)
        norm = plt.Normalize(0, max_runtime)
        bin_colors = cmap(norm(bin_centers))

        bars = ax.bar(bin_centers, hist, width=np.diff(bins), align='center', color=bin_colors, zorder=2, edgecolor='white', linewidth=0.7)

        ax.set_xlabel('Runtime (seconds)', fontsize=12, labelpad=10)
        ax.set_ylabel('Frequency', fontsize=12, labelpad=10)
        ax.set_title('Distribution of Runtimes', fontsize=14, fontweight='bold', pad=15)
        ax.set_xlim(0, max_runtime * 1.02)
        ax.set_ylim(0, max(hist) * 1.1)
        ax.grid(axis='y', alpha=0.75, zorder=0)
        ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.xaxis.set_major_locator(mticker.AutoLocator())
        ax.tick_params(axis='both', labelsize=10, labelcolor = '#333333')


        ax.set_xticks(bins)
        ax.set_xticklabels([f'{b:.2f}' for b in bins], rotation=45, ha='right', fontsize=10)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, label='Runtime (seconds)', orientation='vertical', ax=ax, pad = 0.02, shrink=0.8)
        cbar.ax.tick_params(labelsize=10)
        cbar.outline.set_visible(False)



        plt.tight_layout()
        plt.show()

    else:
        print("No runtime data found.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a histogram of runtimes.")
    parser.add_argument("root_dir", help="The root directory containing the subdirectories with runtime data.")
    args = parser.parse_args()
    plot_runtime_histogram(args.root_dir)
