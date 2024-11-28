import os
import matplotlib.pyplot as plt
import argparse
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

def plot_data_from_directories(root_dir):
    x_labels = []
    y_values = []

    subdir_names = sorted(os.listdir(root_dir))
    for subdir_name in subdir_names:
        subdir_path = os.path.join(root_dir, subdir_name)

        if not os.path.isdir(subdir_path):
            continue

        runtime_file_path = os.path.join(subdir_path, 'runtime.txt')

        if os.path.exists(runtime_file_path):
            try:
                with open(runtime_file_path, 'r') as runtime_file:
                    runtime = float(runtime_file.readline().strip())

                x_labels.append(subdir_name)
                y_values.append(runtime)

            except ValueError as e:
                 print(f"Error processing file in directory {subdir_name}: {e}")
                 continue

    if x_labels and y_values:

        x_indices = np.arange(len(x_labels))

        fig, ax = plt.subplots()

        cmap = LinearSegmentedColormap.from_list("GreenYellowRed", ["green", "yellow", "red"], N=256)

        norm = plt.Normalize(min(y_values), max(y_values))

        scatter = ax.scatter(x_indices, y_values, c=y_values, cmap=cmap, norm=norm, marker='o', zorder=2)

        ax.set_xlabel('Subdirectories (in order)')
        ax.set_ylabel('Runtime')
        ax.set_title('Runtimes')
        ax.grid(True, zorder=0)
        ax.set_xticks(x_indices)
        ax.set_xticklabels(['' for _ in x_labels])

        plt.colorbar(scatter, label='Runtime', orientation="vertical")

        plt.show()

    else:
         print("No data found to plot.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot data from subdirectories.")
    parser.add_argument("root_dir", help="The root directory containing the data subdirectories.")
    args = parser.parse_args()

    root_directory = args.root_dir
    plot_data_from_directories(root_directory)
