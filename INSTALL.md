Building an Example or a Demo of CGAL
=====================================

Since Version 5.0, CGAL is a header-only library, hence it is not needed to build it. Usage of CGAL should simply amount to:

``` {.bash}
git clone https://github.com/CGAL/cgal.git /path/to/cgal.git
cd /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DIR=/path/to/cgal.git ../..
make
```

in the case of building some CGAL-provided examples in debug mode.

Note that although CGAL is a header-only library, some parts of it must link to several external libraries, such as GMP, MPFR, etc.

Building a Program Using CGAL
=============================

If you wish to build a program that is not provided with CGAL and does not already have a `CMakeLists.txt`,
you can trigger the creation of a basic `CMakeLists.txt` by calling the script [`cgal_create_cmake_script`](Scripts/scripts/cgal_create_cmake_script)
found in `/path/to/cgal.git/Scripts/scripts/` at the root of your program directory.

``` {.bash}
 git clone https://github.com/CGAL/cgal.git /path/to/cgal.git
 cd /path/to/your/program
 /path/to/cgal.git/Scripts/scripts/cgal_create_cmake_script
 mkdir -p build/debug
 cd build/debug
 cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DIR:PATH=/path/to/cgal.git ../..
 make your_program
```

Since the basic `CMakeLists.txt` created by the script `cgal_create_cmake_script` cannot
guess which part(s) of CGAL you are using, it does not link with any optional third party 
dependency of CGAL. You should look at the documentation of the package(s) that you
are using to learn which dependencies you must add. The `CMakeLists.txt`
of the examples and demos provided with the package(s) that you are using can be used
to complete your basic `CMakeLists.txt`.


Repository Structure
====================

If you have downloaded a source release instead of cloning the Git repository, the files will be organized in a slightly different way, see the [Layout of the CGAL Git Repository](README.md).

Documentation
=============

For more information see the [CGAL manual](https://doc.cgal.org/latest/Manual/general_intro.html).
