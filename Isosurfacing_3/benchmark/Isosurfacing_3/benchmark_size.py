from benchmark_util import *

# KERNEL_SIMPLE_CARTESIAN_DOUBLE
# KERNEL_SIMPLE_CARTESIAN_FLOAT
# KERNEL_CARTESIAN_DOUBLE
# KERNEL_EPIC
kernel = "KERNEL_CARTESIAN_DOUBLE"

# SCENARIO_GRID_SPHERE
# SCENARIO_IMPLICIT_SPHERE
# SCENARIO_IMPLICIT_IWP
# SCENARIO_SKULL_IMAGE
scenario = "SCENARIO_SKULL_IMAGE"

# TAG_SEQUENTIAL
# TAG_PARALLEL
tag = "TAG_PARALLEL"

# ALGO_MARCHING_CUBES
# ALGO_DUAL_CONTOURING
algorithm = "ALGO_DUAL_CONTOURING"
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
