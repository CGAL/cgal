from benchmark_util import *

scenario = "SCENARIO_GRID_SPHERE"
kernel = "KERNEL_SIMPLE_CARTESIAN_FLOAT"
algorithm = "ALGO_DUAL_CONTOURING"
min_threads = 1
max_threads = 8
cells = 1000000

build(scenario, kernel, algorithm)

data = []

for t in range(min_threads, max_threads):
	time = execute(cells)

	data.append([scenario, kernel, algorithm, t, cells, 0])

df = pd.DataFrame(data, comulmns=["scenario", "kernel", "algorithm", "threads", "cells", "time"])

df.to_csv("benchmark_threads.csv")
