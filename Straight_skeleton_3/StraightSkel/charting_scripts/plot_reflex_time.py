import os
import matplotlib.pyplot as plt
import argparse
from matplotlib.colors import LinearSegmentedColormap

def plot_data_from_directories(root_dir):
    x_values = []
    y_values = []

    for subdir_name in os.listdir(root_dir):
        subdir_path = os.path.join(root_dir, subdir_name)

        if not os.path.isdir(subdir_path):
            continue

        runtime_file_path = os.path.join(subdir_path, 'runtime.txt')
        reflex_edge_count_file_path = os.path.join(subdir_path, 'reflex_edge_count.txt')

        if os.path.exists(runtime_file_path) and os.path.exists(reflex_edge_count_file_path):
            try:
                with open(runtime_file_path, 'r') as runtime_file:
                    runtime = float(runtime_file.readline().strip())

                with open(reflex_edge_count_file_path, 'r') as reflex_edge_count_file:
                    reflex_edge_count = int(reflex_edge_count_file.readline().strip())

                x_values.append(reflex_edge_count)
                y_values.append(runtime)

            except ValueError as e:
                 print(f"Error processing files in directory {subdir_name}: {e}")
                 continue

    if x_values and y_values:
         fig, ax = plt.subplots(figsize=(8, 6))

         cmap = LinearSegmentedColormap.from_list("GreenYellowRed", ["green", "yellow", "red"], N=256)
         norm = plt.Normalize(min(y_values), max(y_values))

         scatter = ax.scatter(x_values, y_values, c=y_values, cmap=cmap, norm=norm, marker='o', zorder = 2)

         ax.set_xlabel('Reflex Edge Count')
         ax.set_ylabel('Runtime (seconds)')
         ax.set_title(f'Runtime vs. Reflex Edge Count (Data Points: {len(x_values)})')  # Add the number of data points to the title
         ax.grid(True, zorder=0)

         # Create a new set of axes for the colorbar legend
         cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
         plt.colorbar(scatter, cax=cbar_ax, label='Runtime')

         # Create a new set of axes for the zoomed-in plot
         zoom_ax = ax.inset_axes([0.55, 0.15, 0.35, 0.35])  # [x0, y0, width, height] in axes coordinates

         # Filter the data for X values below 100
         zoom_x_values = [x for x in x_values if x < 100]
         zoom_y_values = [y for x, y in zip(x_values, y_values) if x < 100]

         # Plot the zoomed-in data
         zoom_scatter = zoom_ax.scatter(zoom_x_values, zoom_y_values, c=zoom_y_values, cmap=cmap, norm=norm, marker='o', zorder=2)

         zoom_ax.set_xlabel('Reflex Edge Count')
         zoom_ax.set_ylabel('Runtime (seconds)')
         zoom_ax.set_title(f'Zoomed-in Plot (X < 100) (Data Points: {len(zoom_x_values)})')  # Add the number of data points to the subplot title
         zoom_ax.grid(True, zorder=0)

         plt.tight_layout(rect=[0, 0, 0.9, 1])  # Adjust the spacing and leave room for the colorbar legend
         plt.show()
    else:
         print("No data found to plot.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot data from subdirectories.")
    parser.add_argument("root_dir", help="The root directory containing the data subdirectories.")
    args = parser.parse_args()

    root_directory = args.root_dir
    plot_data_from_directories(root_directory)
