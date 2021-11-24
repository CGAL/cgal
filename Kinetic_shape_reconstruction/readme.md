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