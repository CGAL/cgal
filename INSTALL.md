Building a Program Using CGAL
=============================

Since Version 5.0, CGAL is header-only, hence it does not create a library. Some parts of it, however, link to several external libraries, such as GMP, MPFR, etc. 

If you are trying to build examples or tests that do not already have a `CMakeLists.txt`,
you can trigger its creation by calling the script [`cgal_create_cmake_script`](Scripts/scripts/cgal_create_cmake_script)
found in `/path/to/cgal.git/Scripts/scripts/` at the root of the example/test directory.
Here is a recipe for the examples of the 2D Triangulation package:

``` {.bash}
 git clone https://github.com/CGAL/cgal.git /path/to/cgal.git
 cd /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
 /path/to/cgal.git/Scripts/scripts/cgal_create_cmake_script
 cd /path/to/cgal.git/Triangulation_2/examples/Triangulation_2
 mkdir -p build/debug
 cd build/debug
 cmake -DCGAL_DIR:PATH=/path/to/cgal.git ../..
 make
```

If, instead of the git repository you downloaded a source release, the files will be organized in a slightly different way, see the [Layout of the CGAL Git Repository](README.md).

For more information see the [CGAL manual](https://doc.cgal.org/latest/Manual/general_intro.html).
