#
# Small helper script to create query files using the classcial Fréchet algorithm
#

import subprocess

dimensions = [3, 100]
eps = 10e-10 # note that this is just set the same as in Compute_classical_Frechet_distance_3.cpp

def read_curve_filenames(d):
    # read file names from dataset
    file = open(f'../curves/generated_{d}d/dataset.txt')

    # read all curves
    filenames = []
    for line in file:
        filenames.append(f'../curves/generated_{d}d/{line.strip()}')

    return filenames

def write_queries(queries, d):
    # open file
    file = open(f'../queries/generated_{d}d.txt', 'w')

    # write queries
    for q in queries:
        # larger query
        file.write(f'{q[0]} {q[1]} {q[2]+2*eps} 1\n')
        # smaller query
        file.write(f'{q[1]} {q[0]} {q[2]-2*eps} 0\n')

def frechet_distance(curve1, curve2):
    args = ('build/Compute_classical_Frechet_distance_{d}', curve1, curve2)
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    return float(output)

for d in dimensions:
    # read curves
    curve_filenames = read_curve_filenames(d)

    # compute pairwise Fréchet distance
    queries = []
    for i in range(len(curve_filenames)):
        for j in range(i+1,len(curve_filenames)):
            dist = frechet_distance(curve_filenames[i], curve_filenames[j])
            queries.append((i,j,dist))

    # write queries to file
    write_queries(queries, d)
