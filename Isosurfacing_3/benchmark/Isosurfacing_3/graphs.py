
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


def add_triangle_graph(data, label, factor):
    x = data["cells"]
    y = data["polygons"] * factor
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


data = pd.read_csv("implicit_iwp_mc_1.csv")
add_threads_graph(data, "MC")

#data = pd.read_csv("threads_grid.csv")
#add_threads_graph(data, "grid")

xt = np.arange(0, min(max(data["threads"]), 9) + 0.1, 2)
if max(data["threads"]) > 10:
    print("more")
    xt = np.concatenate((xt, np.arange(10, max(data["threads"]) + 0.1, 2)))

yt = np.arange(0, 20 + 0.1, 2)

print(xt)
plt.xticks(xt)
plt.yticks(yt)
plot_graph("perf_threads.svg", "", False, "performance [10^3 cubes/s]", "cores")


data = pd.read_csv("size_iwp_mc.csv")
add_size_graph(data, "MC")

data = pd.read_csv("size_iwp_dc.csv")
add_size_graph(data, "DC")

plot_graph("perf_size.svg", "", False, "performance [10^3 cubes/s]", "cells")

data = pd.read_csv("size_iwp_mc.csv")
add_triangle_graph(data, "MC", 1)

data = pd.read_csv("size_iwp_dc.csv")
add_triangle_graph(data, "DC", 1)

plot_graph("triangles_size.svg", "", False, "triangles", "cells")
