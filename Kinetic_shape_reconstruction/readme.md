** GENERAL **

All the test data are placed in the folder examples/data. The name of the folder indicates the type of data or its complexity.
The complexity 1 is low while the complexity 6 is high.

Examples are the entry point to the package:
- kinetic_2d_example: 2D kinetic algorithm on random set of segment
- kinetic_precomputed_shapes_example: 3D kinetic algorithm on a set of user-defined polygons
- kinetic_random_shapes_example: 3D kinetic algorithm on a set of random polygons
- kinetic_reconstruction_example: given a ply with classified point cloud, it first fits polygons using Region Growing, then runs 3D kinetic, and finally filters out exterior volumes creating a reconstructed model

Tests:
- kinetic_2d_stress_test: tests the 2D kinetic algorithm
- kinetic_3d_test_all: tests 3D kinetic algorithm on all data from the examples folder


** INTERSTING CASES **

- EPICK versus EPECK:
By running:
./kinetic_precomputed_shapes_example data/stress-test-4/test-9-rnd-polygons-12-4.off
first with EPICK and then with EPECK shows a huge difference in runtime already. This
data set contains 12 random polygons, each having 4 vertices. The time for EPICK is
milliseconds while for EPECK about 3 minutes. It can also be noticed that most of
the time is spent at the last iterations while the first iterations are very fast.
It is even slower for `Simple_cartesian<Gmpq>`.

It is probably due to the fact that at the latest iterations the computation involves
all operations carried out from the first iteration, which is slow.

Possible solutions:
- Use EPICK and compute always with respect to the input. E.g. intersections between lines should be
done not between lines at the current and previous iterations but between lines at
the current and first iterations. The reason for this optimization is that if we use simply EPICK,
when the number of input polygons grow, we bump into an issue of the accumulating error that gets
bigger and bigger with each iteration and at some point can break the results. It does not happen
with EPECK but we lose speed.
- Use EPECK only when it is absolutely necessary like when computing intersections and
use EPICK for all other computations. This way we avoid accumulating errors and should
keep speed.


** INTERNALS **

- Epsilon usage inside the code.
- File descriptions.
- Important parts of the code.
- Possible failures and warnings.
- Parameters.
