from benchmark_util import *

scenario = "SCENARIO_IMPLICIT_IWP"
kernel = "KERNEL_SIMPLE_CARTESIAN_DOUBLE"
algorithm = "ALGO_MARCHING_CUBES"
tag = "TAG_SEQUENTIAL"
threads = 1
exponent = 1.2
min_cells = 100000
max_cells = 100000000

build(scenario, kernel, algorithm, tag)

data = []

c = min_cells
while c < max_cells:
	n = int(c ** (1.0 / 3.0))

	res = execute(n, threads, 5)
	data.append([scenario, kernel, algorithm, tag, threads, int(c), res["time"], res["polygons"], res["points"]])

	c *= exponent

df = pd.DataFrame(data, columns=["scenario", "kernel", "algorithm", "tag", "threads", "cells", "time", "polygons", "points"])

df.to_csv("benchmark_size.csv")
