from benchmark_util import *

scenario = "SCENARIO_GRID_SPHERE"
kernel = "KERNEL_SIMPLE_CARTESIAN_FLOAT"
algorithm = "ALGO_DUAL_CONTOURING"
threads = 1
exponent = 2
min_cells = 10
max_cells = 1000

build(scenario, kernel, algorithm)

data = []

c = min_cells
while c < max_cells:
	n = int(c ** (1.0 / 3.0))

	time = execute(n)
	data.append([scenario, kernel, algorithm, threads, c, time])

	c *= exponent

df = pd.DataFrame(data, comulmns=["scenario", "kernel", "algorithm", "threads", "cells", "time"])

df.to_csv("benchmark_size.csv")
