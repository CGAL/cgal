Building CGAL Libraries From a Branch
=====================================

Building CGAL using the *branch build* presented here keeps the
build-sources attached to the Git repository.

Note that we do not document here what are the dependancies and cmake options that
are needed to configure CGAL as they are already documented in the
[installation manual](http://doc.cgal.org/latest/Manual/installation.html).

Branch Build of CGAL
====================
The cmake script at the root of the repository is the one to use to
build the CGAL library from a branch. It will collect the list of packages
of the branch and will append their include folder to the include path.
This is main noticable difference with a build using a regular *flat* release.

Here is an example of how to build the library in Debug:
``` {.bash}
git clone https://github.com/CGAL/cgal.git /path/to/cgal.git
cd /path/to/cgal.git
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug ../..
make
```

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
of where you built the library (environment or cmake variable).

Here is an example of how to build in debug the examples from the 3D Triangulations package:

``` {.bash}
 cmake -DCGAL_DIR:PATH=/path/to/cgal.git/build/debug /path/to/cgal.git/Triangulation_3/examples/Triangulation_3
 make
```

If you're trying to build examples or tests that does not already have a `CMakeLists.txt`, you can trigger its creation by calling the script [`cgal_create_cmake_script`](Scripts/scripts/cgal_create_cmake_script) found in `/path/to/cgal.git/Scripts/scripts/` at the root of the example/test directory. Here is an example for the examples of the 2D Triangulation package:

``` {.bash}
 cd /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
 /path/to/cgal.git/Scripts/scripts/cgal_create_cmake_script
 cd -
 cmake -DCGAL_DIR:PATH=/path/to/cgal.git/build/debug /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
 make
```

Note If You Switch Between Branches
===================================
A build may be outdated after an include/dir has been deleted,
switched or even updated. This might lead to compile problems (link
with outdated version). Thus, it is recommended to build CGAL after
each update, switch, merge of a branch (in particular if directories
have been added/deleted, or cpp files have been added, deleted or
altered).


