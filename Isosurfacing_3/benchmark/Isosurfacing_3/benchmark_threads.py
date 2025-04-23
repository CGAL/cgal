from benchmark_util import *

scenario = "SCENARIO_IMPLICIT_IWP"
kernel = "KERNEL_SIMPLE_CARTESIAN_DOUBLE"
algorithm = "ALGO_MARCHING_CUBES"
tag = "TAG_PARALLEL"
min_threads = 1
max_threads = 26
n = 500
cells = n ** 3

build(scenario, kernel, algorithm, tag)

data = []

for t in range(min_threads, max_threads + 1):
	res = execute(n, t, times=5)

	data.append([scenario, kernel, algorithm, tag, t, cells, res["time"], res["bandwidth"]])

df = pd.DataFrame(data, columns=["scenario", "kernel", "algorithm", "tag", "threads", "cells", "time", "bandwidth"])

df.to_csv("benchmark_threads.csv")
