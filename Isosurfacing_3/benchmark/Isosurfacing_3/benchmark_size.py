from benchmark_util import *

scenario = "SCENARIO_GRID_SPHERE"
kernel = "KERNEL_SIMPLE_CARTESIAN_FLOAT"
algorithm = "ALGO_DUAL_CONTOURING"
tag = "TAG_SEQUENTIAL"
threads = 1
exponent = 1.2
min_cells = 1000000
max_cells = 100000000

build(scenario, kernel, algorithm, tag)

data = []

c = min_cells
while c < max_cells:
	n = int(c ** (1.0 / 3.0))

	time = execute(n)
	data.append([scenario, kernel, algorithm, tag, threads, int(c), time])

	c *= exponent

df = pd.DataFrame(data, columns=["scenario", "kernel", "algorithm", "tag", "threads", "cells", "time"])

df.to_csv("benchmark_size.csv")
