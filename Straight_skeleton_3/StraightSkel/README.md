 S t r a i g h t S k e l
=========================

Implementation of Straight Skeleton and Mesh Offsetting in 3-dimensional space.

 Requirements
--------------

* C++ Compiler
  http://gcc.gnu.org/
* CMake build system
  http://www.cmake.org/
* Boost C++ Libraries
  http://www.boost.org/
* Computational Geometry Algorithms Library (CGAL)
  http://www.cgal.org/

Optional:
* SQLite
  http://www.sqlite.org/
* OpenGL Utility Toolkit (GLUT)
  http://freeglut.sourceforge.net/
* Doxygen
  http://www.stack.nl/~dimitri/doxygen/


 Building
----------

CMake is used as build system.
Therefore you may use "ccmake" or "cmake-gui" to configure the build process (build type,
CGAL and third-party libraries directories, etc.), as usual:

```
mkdir build
cd build
cmake-gui ..
make
```

 Running
---------

The example `src/offset_mesh.cpp` is the main entry point: it is documented pipeline
that checks the sanity of the input and the constructed output, and calls the main
function `CGAL::face_offset()`.

A typical call would be:

```sh
./offset_mesh ~/Data/Customers/XXXX/offset_data_july_2024/input_2105666_0.ply --weight-path ~/Data/Customers/XXXX/offset_data_july_2024/offsets_2105666.txt --save-path last_run/2105666_0 --save-offsets 1
```

See the file for documentation of its parameters, and the documentation of `CGAL::face_offset()`.

Some scripts can be used to automatically run tests:
- `test.sh`, to run a few dozen unit tests
- `run_and_compare.sh`, which takes a path to a directory with triplets of file of the format:
  - "input_XXX.ply" the input
  - "offsets_XXX.txt" the weights along cardinal direction, see the function `assign_weights()` in `offset_mesh.cpp`
  - "output_XXX.ply" an output to compare to.
  and construct the offset mesh and performs some comparisons with the provided result.

 History
---------

StraightSkel was started by Gernot Walzl in the years 2011, 2012, 2013 (MIT).

3D straight skeletons and mesh offsetting were further developed by GeometryFactory
since 2024 (GPL-3.0-or-later OR LicenseRef-Commercial).
