NOTICE
======

Since Version 5.0, CGAL is a header-only library it is not needed
to build and install it. Usage of CGAL should thus simply amount to:

``` {.bash}
git clone https://github.com/CGAL/cgal.git /path/to/cgal.git
cd /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DIR=/path/to/cgal.git
make
```

in the case of the building of an example in debug mode.

For more information head over to the [CGAL manual](https://doc.cgal.org/latest/Manual/general_intro.html).
Note that this page describes the setting of CGAL as a sources release and, as such,
files are organized in a slightly different way, see the [Layout of the CGAL Git Repository](README.md).


Building a Program Using CGAL
=============================

To compile a program using CGAL, simply set `CGAL_DIR` to the location
of the directory containing `CGALConfig.cmake` (for example the root 
of the extracted source archive or the root of a git checkout).

Here is an example of how to build in debug the examples from the 3D Triangulations package:

``` {.bash}
 cd /path/to/cgal.git/Triangulation_3/examples/Triangulation_3
 mkdir -p build/debug
 cd build/debug
 cmake -DCGAL_DIR:PATH=/path/to/cgal.git ../..
 make
```

If you are trying to build examples or tests that do not already have a `CMakeLists.txt`,
you can trigger its creation by calling the script [`cgal_create_cmake_script`](Scripts/scripts/cgal_create_cmake_script)
found in `/path/to/cgal.git/Scripts/scripts/` at the root of the example/test directory.
Here is an example for the examples of the 2D Triangulation package:

``` {.bash}
 cd /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
 /path/to/cgal.git/Scripts/scripts/cgal_create_cmake_script
 cd /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
 mkdir -p build/debug
 cd build/debug
 cmake -DCGAL_DIR:PATH=/path/to/cgal.git ../..
 make
```

For more information head over to the [CGAL manual](https://doc.cgal.org/latest/Manual/general_intro.html).
