NOTICE
======

Since Version 5.0, CGAL is now header-only by default, meaning that you do not need to build and install CGAL. Usage of CGAL as a header-only library
simply amounts to, for example:

``` {.bash}
cd /path/to/cgal/examples/Triangulation_2
mkdir -p build/debug
cd build/debug
cmake -DCMAKE_BUILD_TYPE=Debug -DCGAL_DIR=/path/to/cgal
make
```

For more information head over to the [CGAL manual](https://doc.cgal.org/latest/Manual/general_intro.html).
