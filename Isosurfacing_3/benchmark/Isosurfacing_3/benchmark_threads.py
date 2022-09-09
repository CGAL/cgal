from benchmark_util import *

scenario = "SCENARIO_IMPLICIT_SPHERE"
kernel = "KERNEL_SIMPLE_CARTESIAN_DOUBLE"
algorithm = "ALGO_MARCHING_CUBES"
tag = "TAG_PARALLEL"
min_threads = 1
max_threads = 12
cells = 300

build(scenario, kernel, algorithm, tag)

data = []

for t in range(min_threads, max_threads):
	time = execute(cells, t)

	data.append([scenario, kernel, algorithm, tag, t, cells, time])

df = pd.DataFrame(data, columns=["scenario", "kernel", "algorithm", "tag", "threads", "cells", "time"])

df.to_csv("benchmark_threads.csv")
