 S t r a i g h t S k e l
=========================

StraightSkel is an implementation of Straight Skeleton in 2- and 3-dimensional space.
It is used to compute offsets of closed polygons and polyhedrons.

The straight skeleton is always interior to the shape, and the offsetting is always inwards.

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
Therefore you may use "ccmake" or "cmake-gui" to configure the build process, as usual:

```
mkdir build
cd build
cmake .. # + setup CGAL dir, build type, third party libs, etc.
make
```

 Running
---------

For outward offset mesh generation, a  function is provided, which handles creating
an outside bounding box, inverting the mesh, etc. See sample code: src/offset_mesh.cpp

 Development
-------------

To avoid segmentation fault errors, standard C pointers are used very rarely.
A slightly modified version of the STL's Smart Pointer may be used instead.
(Modification: Print the stack trace in case of an error.)

Organization of Source Code:

src/
  algo/     Algorithms
  data/     Data Structures
  db/       Database
  kernel/   Geometric Kernel (double precision)
  ui/       User Interface
    gl/     OpenGL
    ps/     PostScript
  util/     Tools
  misc/     Miscellaneous (pre- and post-processing helpers for mesh offsetting)
test/       Unit Test Cases

 History
---------

StraightSkel was started by Gernot Walzl in the years 2011, 2012, 2013 (MIT).

StraightSkel was further developed by GeometryFactory since 2024 (GPL-3.0-or-later OR LicenseRef-Commercial).
