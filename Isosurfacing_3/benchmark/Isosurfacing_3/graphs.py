
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def save_svg(file):
    fig = plt.gcf()
    fig.set_size_inches((10, 5), forward=False)
    plt.savefig(file, bbox_inches="tight")


def add_threads_graph(data, label):
    x = data["threads"]
    y = data["cells"] / data["time"] / 10 ** 3
    plt.plot(x, y, label=label)
    plt.legend()


def add_size_graph(data, label):
    x = data["cells"]
    y = data["cells"] / data["time"] / 10 ** 3
    plt.plot(x, y, label=label)
    plt.legend()


def plot_graph(file, name, log, ylabel, xlabel):
    plt.title(name)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if log:
        plt.xscale("log")
    plt.gca().yaxis.grid(color='#cccccc')
    plt.gca().xaxis.grid(color='#cccccc')
    plt.ylim(ymin=0)
    save_svg(file)
    plt.show()


latex_export = True
if latex_export:
    #plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams['font.size'] = "17"


data = pd.read_csv("threads_implicit.csv")
add_threads_graph(data, "implicit")

data = pd.read_csv("threads_grid.csv")
add_threads_graph(data, "grid")

plt.xticks(np.arange(1, max(data["threads"]) + 0.1, 1))
plot_graph("perf_threads.svg", "", False, "performance [10^3 cubes/s]", "cores")


data = pd.read_csv("size_implicit.csv")
add_size_graph(data, "implicit")

data = pd.read_csv("size_grid.csv")
add_size_graph(data, "grid")

plot_graph("perf_size.svg", "", True, "performance [10^3 cubes/s]", "cubes")
