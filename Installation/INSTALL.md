INTRODUCTION
============

This file describes how to install CGAL. The instructions in this file
are for the most common use cases, and cover the command line tools.

For further information, or in case of problems, please see the
detailed installation instructions, which can be found in this
distribution in the file ./doc_html/index.html or on the CGAL website
http://doc.cgal.org/latest/Manual/installation.html

The documentation of CGAL is available in PDF and HTML formats.
It is not bundled with the software but can be downloaded separately
at <http://www.cgal.org/Manual>.

For more information about CGAL, see the <http://www.cgal.org/>.

In the current file, x.y is an implicit replacement for the current version
of CGAL (3.5.1, 3.6, and so on).


PREREQUISITES
=============

To install CGAL, you need 'cmake' and several third-party libraries.
Some are essential for entire CGAL, some are mandatory for particular
CGAL packages, some are only needed for demos.

   * CMake (>= 2.8.11), the build system used by CGAL
     Required for building CGAL

   * Boost (>= 1.48)
     Required for building CGAL and for applications using CGAL
     Required compiled Boost library: Boost.Thread, Boost.System
     Optional compiled Boost library: Boost.Program_options
     http://www.boost.org/   or   http://www.boostpro.com/products/free/
     You need the former if you plan to compile the boost libraries yourself,
     for example because you target 64 bit applications for XP64

   * Exact Arithmetic
     CGAL combines floating point arithmetic with exact arithmetic, in order
     to be fast and reliable. CGAL offers support for GMP and MPFR, for LEDA
     exact number types, as well as a built-in exact number type used when
     none of the other two is installed.
     Required by several examples which have hard coded the number type.

     - GMP (>= 4.1.4)
       http://gmplib.org/
       or precompiled version that can be downloaded with CGAL-x.y-Setup.exe
       based on http://fp.gladman.plus.com/computing/gmp4win.htm

     - MPFR (>= 2.2.1)
       http://www.mpfr.org/
       or precompiled version that can be downloaded with CGAL-x.y-Setup.exe
       based on http://fp.gladman.plus.com/computing/gmp4win.htm

     - LEDA (>= 6.2)
       http://www.algorithmic-solutions.com/leda/index.htm

   * Visualization
     Required for most demos

     - Qt3 (>= 3.3)
       ftp://ftp.qt.nokia.com/qt/source/

     - Qt5 (>= 5.3)
       http://qt-project.org/

     - libQGLViewer
       http://www.libqglviewer.com/

     - Geomview
       http://www.geomview.org/
       Not supported with Visual C++

   * Numerical Libraries
     - EIGEN (>=3.1)
       Required by the packages:
       o Estimation of Local Differential Properties of Point-Sampled Surfaces
       o Approximation of Ridges and Umbilics on Triangulated Surface Meshes
       o Planar Parameterization of Triangulated Surface Meshes
       o Surface Reconstruction from Point Sets
       http://eigen.tuxfamily.org

     - BLAS, LAPACK, ATLAS
       Required by the packages (if EIGEN is not available):
       o Estimation of Local Differential Properties of Point-Sampled Surfaces
       o Approximation of Ridges and Umbilics on Triangulated Surface Meshes
       o Planar Parameterization of Triangulated Surface Meshes
       http://www.netlib.org/blas/, http://www.netlib.org/lapack/
       or precompiled version that can be downloaded with CGAL-x.y-Setup.exe

     - MPFI
       Required by the package:
       o Algebraic Kernel
       https://gforge.inria.fr/projects/mpfi/
       (or shipped with RS http://vegas.loria.fr/rs/)

     - RS (root isolation)
       Required by the package:
       o Algebraic Kernel
       http://vegas.loria.fr/rs/

     - NTL (Number Type Theory)
       Optional for the packages:
       o Polynomial
       o Algebraic Kernel
       http://www.shoup.net/ntl/

   * Miscellaneous

     - zlib
       Optional for the package:
       o Surface Mesh Generator can read compressed images directly
       http://www.zlib.net/

     - ESBTL
       Optional to read PDB files:
       o Import data from a PDB file as CGAL points or weighted points.
         An example is given in package Skin_surface (see example skin_surface_pdb_reader.cpp)
       http://esbtl.sourceforge.net/

CONFIGURATION
=============

To configure CGAL, type
```
  cmake .
```
in the directory that contains this INSTALL file. You can add several options
to this command. The most important ones are

* `-DCMAKE_INSTALL_PREFIX=<dir>`          installation directory [/usr/local]
* `-DCMAKE_BUILD_TYPE=<Debug|Release>`    build type [Release]
* `-DBUILD_SHARED_LIBS=<TRUE|FALSE>`      shared or static libraries [TRUE]
* `-DCMAKE_C_COMPILER=<program>`          C compiler [gcc]
* `-DCMAKE_CXX_COMPILER=<program>`        C++ compiler [g++]

In case you want to add additional compiler and linker flags, you can use

* `-DCGAL_CXX_FLAGS`                      additional compiler flags
* `-DCGAL_MODULE_LINKER_FLAGS`            add. linker flags (static libraries)
* `-DCGAL_SHARED_LINKER_FLAGS`            add. linker flags (shared libraries)
* `-DCGAL_EXE_LINKER_FLAGS`               add. linker flags (executables)

Variants with the additional suffix "_DEBUG" and "_RELEASE" allow to set
separate values for debug and release builds. In case you do not want to add
additional flags, but to override the default flags, replace "CGAL" by
"CMAKE" in the variable names above.

By default demos and examples are not configured. If you want to configure
them at the same time as the CGAL library, you can use

*  `-DWITH_examples=true`
*  `-DWITH_demos=true`

Note that CMake maintains a cache name `CMakeCache.txt`. If you change options
(or your environment changes), it is best to remove that file to avoid
problems.


BUILDING
========

To build the CGAL libraries, type
```
  make
```
(or nmake in a Windows command prompt).
If you want, you can install the CGAL header and libraries. To do so, type
```
  make install
```
You can build all demos or examples by typing
```
  make demos
  make examples
```
If you are interested in the demos or examples of just a particular module,
you can build them in the following way:
```
  make -C demo/Alpha_shapes_2        (or: cd demo/Alpha_shapes_2; make)
  make -C examples/Alpha_shapes_2    (or: cd examples/Alpha_shapes_2; make)
```
A list of all available make targets can be obtained by
```
  make help
```

OUT-OF-SOURCE BUILDS
====================

The above instructions build the CGAL library in the same directory tree as
the CGAL sources. Sometimes it is advisable to place all the generated files
somewhere else. For example, if you want to build the library in several
configurations (debug and release, different compilers, and so on). Using
different build directories keeps all the generated files separated for each
configuration.

In the following, `$CGAL_SRC` denotes the directory with the CGAL sources;
`$CGAL_BUILD` is an arbitrary directory where the generated files will be
placed. You can perform an out-of-source build as follows:
```
  mkdir $CGAL_BUILD
  cd $CGAL_BUILD
  cmake [options] $CGAL_SRC
  make
  make install                       (if desired)
  make demos                         (if desired)
  make examples                      (if desired)
```
Basically, the only difference is the last parameter of the `cmake` command:
`$CGAL_SRC` instead of `.` .

