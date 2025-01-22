#
# Small helper script to create goal directed curves in arbitrary dimensions.
#

import random
import numpy as np
import numpy.linalg as la

import matplotlib
import matplotlib.pyplot as plt

# dimensionality of the curves
D = 2
# approximate average distance between two consecutive vertices
verbose = False
plot = True # plots only the first two dimensions
# number of curves and how many different start and end points
num_curves = 6
num_start_end = 3
start_mean = np.array(D*[0.])
end_mean = np.array(D*[50.])
# random seed
seed = 1237
# curve directory
curve_name_prefix = f'../curves/generated_{D}d/'

rng = np.random.default_rng(seed)

def norm(x):
    return la.norm(x, np.inf)

def create_random_point(mean=np.zeros(D), stddev=1.3):
    return rng.normal(mean, stddev)

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
        start = rng.choice(starts)
        end = rng.choice(ends)
        curves.append(create_curve(start, end))
    return curves

def write_to_file(curves, prefix = "curve"):
    for i in range(len(curves)):
        f = open(prefix + str(i) + ".txt", "w")
        for p in curves[i]:
            if D > 3:
                f.write(f'{D}') # only Point_d has the dimension of the point at the beginning
            for i in range(len(p)):
                if i > 0 or D > 3:
                    f.write(" ")
                f.write(str(p[i]))
            f.write("\n")

starts = []
ends = []
for i in range(num_start_end):
    starts.append(create_random_point(start_mean, 8.))
    ends.append(create_random_point(end_mean, 8.))

curves = create_curves(starts, ends, num_curves)

# print
if verbose:
    print(curves)

# plot
if plot:
    for i, curve in enumerate(curves):
        X = [p[0] for p in curve]
        Y = [p[1] for p in curve]
        if i == 5:
            plt.plot(X, Y, linewidth=5)
        elif i == 1 or i == 3:
            plt.plot(X, Y, linewidth=2)
        elif i == 0:
            plt.plot(X, Y, color="gray", linewidth=2)
        else:
            plt.plot(X, Y, color="gray", alpha=.3)
    plt.axis('scaled')
    plt.tight_layout()
    plt.show()

# export
write_to_file(curves, curve_name_prefix)
