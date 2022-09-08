from benchmark_util import *

scenario = "SCENARIO_GRID_SPHERE"
kernel = "KERNEL_SIMPLE_CARTESIAN_FLOAT"
algorithm = "ALGO_MARCHING_CUBES"
tag = "TAG_PARALLEL"
min_threads = 1
max_threads = 8
cells = 500

build(scenario, kernel, algorithm, tag)

data = []

for t in range(min_threads, max_threads):
	time = execute(cells)

	data.append([scenario, kernel, algorithm, tag, t, cells, time])

df = pd.DataFrame(data, columns=["scenario", "kernel", "algorithm", "tag", "threads", "cells", "time"])

df.to_csv("benchmark_threads.csv")
