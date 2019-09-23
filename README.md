[![Build Status](https://travis-ci.org/CGAL/cgal.svg?branch=master)](https://travis-ci.org/CGAL/cgal)

![CGAL](Installation/doc_html/images/cgal_2013_grey.png)

The Computational Geometry Algorithms Library (CGAL) is a C++ library that
aims to provide easy access to efficient and reliable algorithms in
computational geometry.

CGAL Releases
=============
The primary vector of distribution of CGAL are sources tarballs, released
twice a year, announced on [the web site of CGAL](https://www.cgal.org/).

Getting Started with CGAL
=========================

**Since version 5.0, CGAL is header-only by default, meaning that
it is no longer necessary to build or install CGAL before it can be used.**

Head over to the [CGAL manual](https://doc.cgal.org/latest/Manual/general_intro.html)
for usage guides and tutorials that will get you started smoothly.

Compilation and Installation
============================

If you do not wish to use the header-only mode of CGAL, it is still possible to build CGAL itself.

The compilation and installation of CGAL from a sources tarball are
described in the [CGAL installation manual](https://doc.cgal.org/latest/Manual/installation.html)
and in the file [INSTALL.md](Installation/INSTALL.md) that is at the root
of any sources tarball.

You can also clone this repository and compile from your local Git repository.
This kind of compilation is called a *branch build*, and is
described in the file [INSTALL.md](INSTALL.md) that is at the root of the
Git repository.

License
=======
See the file [LICENSE.md](LICENSE.md).

CGAL Git Repository Layout
==========================

The Git repository of CGAL has a different layout from release tarballs. It
contains a `CMakeLists.txt` file that serves as anchor for configuring and building programs,
and a set of subfolders, so called *packages*. Most packages
implement a data structure or an algorithm for CGAL (e.g., `Convex_hull_2`,
or `Triangulation_3`); however some packages serve special needs:

* `Installation` - meta-files and CMake-support
* `Maintenance` - infrastructural support
* `Core`, `CGALimageIO`, `Qt_widget`, `GraphicsView` - component libraries
* `Scripts` - scripts to simplify developer's and user's work
* `Testsuite` - infrastructure for testsuite
* `Documentation` - infrastructure for CGAL's manual
* `STL_Extension` - extensions to the standard template library

More Information
================
* [The CGAL web site](https://www.cgal.org/)
* [Latest CGAL release documentation pages](https://doc.cgal.org/)
* [Latest CGAL master documentation pages, updated once a week](https://cgal.geometryfactory.com/CGAL/doc/master/)
* [CGAL daily testsuite results](https://cgal.geometryfactory.com/CGAL/testsuite/)
* [Guidelines for CGAL developers](https://github.com/CGAL/cgal/wiki/Guidelines) and [Informations for new developers](https://github.com/CGAL/cgal/wiki/Information-for-New-Developers)
