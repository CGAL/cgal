[![Build Status](https://travis-ci.org/CGAL/cgal.svg?branch=master)](https://travis-ci.org/CGAL/cgal)

#![CGAL](Installation/doc_html/images/cgal_2013_grey.png)

The Computational Geometry Algorithms Library (CGAL) is a C++ library that
aims to provide easy access to efficient and reliable algorithms in
computational geometry.

CGAL releases
=============
The primary vector of distribution of CGAL are sources tarballs, released
twice a year, announced on [the web site of CGAL](http://www.cgal.org/).
The sources distributed that way can be built using the
[CGAL installation manual](http://doc.cgal.org/latest/Manual/installation.html).

CGAL Git repository layout
==========================

The Git repository of CGAL has a different layout from release tarballs. It
contains a `CMakeLists.txt` file that serves as anchor for building, and a
set of subfolders, so called *packages*. Most packages
implement a data structure or an algorithm for CGAL (e.g., `Convex_hull_2`,
or `Triangulation_3`); however some packages serve special needs:

* `Installation` - meta-files and CMake-support
* `Maintenance` - infrastructural support
* `Core`, `CGALimageIO`, `Qt_widget`, `GraphicsView` - component libraries
* `Scripts` - scripts to simplify developer's and user's work
* `Testsuite` - infrastructure for testsuite
* `Documentation` - infrastructure for CGAL's manual
* `STL_Extension` - extensions to the standard template library

Compilation and installation
============================
The compilation and installation of CGAL from a sources tarball is
described in the
[CGAL installation manual](http://doc.cgal.org/latest/Manual/installation.html)
and in the file [INSTALL.md](Installation/INSTALL.md) that is at the root
of any sources tarball.

CGAL developers, however, usually compile CGAL directly from a local Git
repository. That kind of compilation is called a *branch build*, and is
described in the file [INSTALL.md](INSTALL.md) that is at the root of the
Git repository.

License
=======
See the file [LICENSE.md](LICENSE.md).

More information
================
* [The CGAL web site](http://www.cgal.org/)
* [Latest CGAL release documentation pages](http://doc.cgal.org/)
* [Latest CGAL master documentation pages, updated once a week](https://cgal.geometryfactory.com/CGAL/doc/master/)
* [CGAL daily testsuite results](https://cgal.geometryfactory.com/CGAL/testsuite/)
