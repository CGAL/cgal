NOTICE
======

Since Version 5.0, CGAL is now header-only by default, meaning that it is no longer mandatory
to build (and install) CGAL and its libraries. Usage of CGAL as a header-only library should
thus simply amount to:

``` {.bash}
git clone https://github.com/CGAL/cgal.git /path/to/cgal.git
cd /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DIR=/path/to/cgal.git
make
```

in the case of the building of an example in debug mode.

More information on dependencies, configuration flags, and useful scripts to build programs that are not shipped
with CGAL can be found on header-only usage at https://doc.cgal.org/latest/Manual/general_intro.html,
noting that this page describes the setting of CGAL as a sources release and, as such,
files are organized in a slightly different way, see the [Layout of the CGAL Git Repository](README.md).

BRANCH BUILD OF CGAL
====================

Although not recommended, it is still possible to build (and install) CGAL.

The cmake script at the root of the repository is the one to use to
build the CGAL library from a branch. It will collect the list of packages
of the branch and will append their include folder to the include path.
This is main noticeable difference with a build using a regular *flat* release.

Here is an example of how to build the library in Release:
``` {.bash}
git clone https://github.com/CGAL/cgal.git /path/to/cgal.git
cd /path/to/cgal.git
mkdir -p build/release
cd build/release
cmake -DCMAKE_BUILD_TYPE=Release ../..
make
```

Note that *no installation is required* to use that version of CGAL once it has been compiled.

Building a Program Using CGAL
=============================

To compile a program using CGAL, simply set `CGAL_DIR` to the location
of where you have built the library (environment or cmake variable).

Here is an example of how to build in debug the examples from the 3D Triangulations package:

``` {.bash}
 cd /path/to/cgal.git/Triangulation_3/examples/Triangulation_3
 mkdir -p build/debug
 cd build/debug
 cmake -DCGAL_DIR:PATH=/path/to/cgal.git/build/debug ../..
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
 cmake -DCGAL_DIR:PATH=/path/to/cgal.git/build/debug ../..
 make
```

More information on building and installing CGAL is available in the CGAL manual,
at https://doc.cgal.org/latest/Manual/installation.html, noting that this page describes the setting
of CGAL as a sources release and, as such, files are organized in a slightly different way,
see the [Layout of the CGAL Git Repository](README.md).

Note If You Switch Between Branches
===================================

A build may be outdated after an include/dir has been deleted,
switched or even updated. This might lead to compile problems (link
with outdated version). Thus, it is recommended to build CGAL after
each update, switch, merge of a branch (in particular if directories
have been added/deleted, or cpp files have been added, deleted or
altered).


