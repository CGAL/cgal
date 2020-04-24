# :heavy_exclamation_mark: Notice :heavy_exclamation_mark:
This repository is only meant for reviewing purposes. Please do not distribute or link to this private repository. The final implementation will be open-source and integrated in the [Computational Geometry Algorithm Library (CGAL)](https://github.com/CGAL/cgal/).

<p align="center"><img src="GPDT2.jpg" alt="drawing" width="500"/></p>

# Delaunay Triangulations of Lattice Points

This repository contains the implementation of the algorithm described in XXXXXXXXXXXXXXXXXX, currently submitted to [ESA 2020](http://algo2020.di.unipi.it/ESA2020/index.html).

# Structure of the Repository

Our implementation enhances and builds on top of existing triangulations of CGAL. Consequently, the repository follows the same structure as the git repository of CGAL, with each component of CGAL being a directory. Our contribution mainly resides within the `Periodic_2_triangulation_2` and `Periodic_3_triangulation_3` packages: a completely new triangulation class is introduced, along with its vertex, face, and cell classes.

The CGAL codebase was forked from commit https://github.com/CGAL/cgal/commit/2167d4ffef72241e382d0779b448835660af6ca9 (this can be used to get a complete set of changes).

# Usage

Two examples illustrate our algorithm, [generic_p2t2.cpp](https://github.com/MaelRL/cgal/blob/Generic_P2T2/Periodic_2_triangulation_2/examples/Periodic_2_triangulation_2/generic_p2t2.cpp) and [generic_p3t3.cpp](https://github.com/MaelRL/cgal/blob/Generic_P2T2/Periodic_3_triangulation_3/examples/Periodic_3_triangulation_3/generic_p3t3.cpp).

Downloading this code, compiling and executing the examples can be achieved as follows:
```
mkdir /path/to/cgal
git clone git@github.com:GPT-Authors/Periodic-Delaunay-Triangulations-on-Lattices.git /path/to/cgal // or download the code directly from https://github.com/GPT-Authors/Periodic-Delaunay-Triangulations-on-Lattices.git
cd /path/to/cgal/Periodic_2_triangulation_2/examples/Periodic_2_triangulation_2
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DIR=/path/to/cgal ../..
make generic_p2t2
```

See also https://github.com/CGAL/cgal/blob/master/INSTALL.md and CGAL's documentation for further installation instructions, if necessary.

A typical output will be:
```
./generic_p2t2 [number_of_vertices]
random seed: 5588160 // the seed used to initialize the random number generator
Basis vectors:
24.7802 -0.222676 -50.4873
-10.9504 -8.12726 -0.595036
-59.6871 -37.184 -41.2713
Inserting 10000 random points...
Transition to Phase 2 after 890 vertices
Done!
Number of vertices: 10000
Number of cells: 67558
Is the triangulation 1-cover? true // whether we transitioned to Phase 2 at some point or stayed in Phase 1
Time: 0.650457 s
```

As well as creating a file called `final.off`, which can be opened with a large number of mesh visualization softwares such as the CGAL Polyhedron Demo, Meshlab, etc.

A specific lattice can easily be specified through a set of 2 (3 in 3D) independent vectors and recompiling the program.

# Functionality

Although the core functionalities (location of a point, insertion of a point) of our new periodic triangulations are present, some additional functionality such as point removal and weighted triangulations are not yet implemented. Those functionalities will be added in the final version before its integration in CGAL.

# License

The code is currently private and has no license. The final license will be GPL, which is the license of the CGAL Periodic Triangulations packages.
