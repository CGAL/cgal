NOTICE
======

Since Version 5.0, CGAL is now header-only by default, meaning that it is no longer mandatory
to build (and install) CGAL and its libraries. Usage of CGAL as a header-only library
simply amounts to, for example:

``` {.bash}
cd /path/to/cgal/examples/Triangulation_2
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DIR=/path/to/cgal
make
```

More information on header-only dependencies, configuration flags, and useful scripts to build programs
that are not shipped with CGAL is available in this distribution in the file ./doc_html/index.html
("Getting Started with CGAL"), or on the CGAL website at https://doc.cgal.org/latest/Manual/general_intro.html

BUILDING AND INSTALLING CGAL
============================

Although not recommended, it is still possible to build and install CGAL. More information on
building and installing CGAL is available in this distribution in the file ./doc_html/index.html
("Advanced" > "Installation"), or on the CGAL website at https://doc.cgal.org/latest/Manual/installation.html

