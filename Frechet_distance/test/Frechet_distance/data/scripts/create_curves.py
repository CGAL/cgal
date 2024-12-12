#
# Small helper script to create goal directed curves in arbitrary dimensions.
#

import random
import numpy as np
import numpy.linalg as la

import matplotlib
import matplotlib.pyplot as plt

# dimensionality of the curves
D = 100
# approximate average distance between two consecutive vertices
verbose = False
plot = False # plots only the first two dimensions
# number of curves and how many different start and end points
num_curves = 10
num_start_end = 3
start_mean = np.array(D*[0.])
end_mean = np.array(D*[20.])
# random seed
seed = 42
# curve directory
curve_name_prefix = "../curves/generated_100d/"

rng = np.random.default_rng(seed)

def norm(x):
    return la.norm(x, np.inf)

def create_random_point(mean=np.zeros(D)):
    return rng.normal(mean, 1.3)

def create_curve(start, end):
    curve = [start]
    while norm(end - curve[-1]) > np.sqrt(D):
        step = curve[-1] + (end - curve[-1])/norm(end - curve[-1])
        noisy_step = create_random_point(step)
        curve.append(noisy_step)
    curve.append(end)
    return curve

def create_curves(starts=None, ends=None, num_curves=1):
    curves = []
    for i in range(num_curves):
        start = random.choice(starts)
        end = random.choice(ends)
        curves.append(create_curve(start, end))
    return curves

def write_to_file(curves, prefix = "curve"):
    for i in range(len(curves)):
        f = open(prefix + str(i) + ".txt", "w")
        for p in curves[i]:
            for i in range(len(p)):
                f.write(str(p[i]))
                if i < len(p)-1:
                    f.write(" ")
            f.write("\n")

starts = []
ends = []
for i in range(num_start_end):
    starts.append(create_random_point(start_mean))
    ends.append(create_random_point(end_mean))

curves = create_curves(starts, ends, num_curves)

# print
if verbose:
    print(curves)

# plot
if plot:
    for curve in curves:
        X = [p[0] for p in curve]
        Y = [p[1] for p in curve]
        plt.plot(X, Y)
    plt.show()

# export
write_to_file(curves, curve_name_prefix)
