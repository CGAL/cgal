Building CGAL Libraries and Executables based on CGAL
=====================================================

Building CGAL using the *branch build* presented here keeps the
build-sources attached to the Git repository. A developer, let's call her
Jenny, can build CGAL from any CGAL branch accessible through a proper
*source code management* (SCM) client, such as Git.

This way Jenny can modify CGAL sources and test them without disturbing
other developers and without disturbance from other developers, and at the
same time keep (modified) CGAL sources in sync with an SCM repository and
in turn with other developers. This documentation mainly discusses how to
generate Unix Makefiles with CMake; however other generators are supported
by CMake, too. In general, it is also advised to be familiar with the usual
build procedure as the branch-build only selects another path to the main
`CMakeLists.txt` file (i.e., picks it "from a branch" instead of "from a
sources tarball").

Example
=======

Assume that Jenny wants to try the new cool feature developed by Adam
about coloring 2D triangulations. Adam is developing his feature in the
branch `Triangulation_2-coloring-adam`. Jenny got a working copy of
`next` (or of a different branch, e.g., any of her feature branches),
e.g., by calling

``` {.bash}
 > cd /path/to/my/cgal_repository/
 > git checkout Triangulation_2-coloring-adam
 > git branch
 master
 * Triangulation_2-coloring-adam
 ...
```

Jenny aims to compile a demo in that branch using the headers from the
branch and the CGAL libraries built from the sources contained in the
branch. Learn next how Jenny does it.

Branch Build of CGAL
====================

Jenny has many options to use the *branch build* to build CGAL and then
compile some examples. The build is categorized either as in-source or
out-source. The former means that generated directories, generated
files, and the sources are placed side by side. The latter means that
generated directories and files are placed under a completely different
directory, which, for example, can be removed when not needed without
any impact on the sources. Among the different options, the in-source
option is, while possible, less typical. Thus, it is not described here.

Jenny can choose among three typical use-cases. Each case has advantages
and disadvantages, and the preferred choice depends on how Jenny's daily
work looks like.

Preamble: Options for `cmake`
-----------------------------

All three have in common that Jenny might enhance the calls to `cmake`
with platform-specific options or other flags as switching on or off
components of CGAL, e.g.,

``` {.bash}
 -DWITH_CGAL_Core=ON -DWITH_CGAL_Qt3=OFF -DWITH_CGAL_Qt4=ON -DWITH_CGAL_ImageIO=ON -DWITH_examples=OFF -DWITH_demos=OFF
```

or configuring external libraries, e.g.,

``` {.bash}
 -DWITH_MPFI=ON -DWITH_RS=ON -DWITH_LEDA=ON -DLEDA_CXX_FLAGS=-ffriend-injection -DLEDA_INCLUDE_DIR=$EXACUS_LEDA/incl/ -DLEDA_LIBRARIES=$EXACUS_LEDA/libleda.so -DLEDA_LINKER_FLAGS=-lX11  
```

or any other option valid for `cmake`, e.g.

``` {.bash}
 -DCMAKE_BUILD_TYPE=Debug
```

All possible options and external libraries can be found in the
[Installation
manual](http://doc.cgal.org/latest/Manual/installation.html).
So, assume next that Jenny knows the options according to her platform
and liking when calling `cmake`. For simplicity, these options are not
mentioned in the next examples.

Those options can also be set using the CMake GUI `cmake-gui`.

Building locally to the working copy
------------------------------------

Jenny can set `$CGAL_DIR` to `/path/to/my/cgal_repository/build`. It
couples the sources and the build closer while still having an
out-of-source-build.

``` {.bash}
  > mkdir /path/to/my/cgal_repository/build
  > cd /path/to/my/cgal_repository/build
  > cmake ..
  > make
  > export CGAL_DIR=/path/to/my/cgal_repository/build
```

Using a single version of CGAL
------------------------------

For many packages there is no trace of the package code in the CGAL
library objects, as the entire code of the packages resides only in
header files. As a matter of fact, the library objects themselves are
very thin and thus building CGAL is a quick process and the generated
library objects consume little space. In many cases it is sufficient to
maintain a single version of CGAL at any time, for example, when only a
single feature is developed, when several features are developed in
parallel but the different corresponding branches can share a common
version of CGAL, or when several features are developed in a single
working directory and the corresponding branches are switched rarely. In
such cases the environment variable <CGAL_DIR>, which points to the
build target directory, is set once and never changes.

Assume that Jenny would like to build CGAL from her local
repositorylocated at `/path/to/my/cgal_repository` and place the result
of the build in `~/CGAL/build`. She issue the following commands:

``` {.bash}
 > export CGAL_DIR=$HOME/CGAL/build
 > cd $CGAL_DIR
 > cmake /path/to/my/cgal_repository
 > make
```

Using multiple versions of CGAL
-------------------------------

In this setup, Jenny develops several features and often switches
between the corresponding branches in repository. She would like to keep
a build version of CGAL for each of the features she is developing.
Thus, for a branch `Triangulation_2-colored-adam` she builds CGAL in
`~/CGAL/builds/Triangulation_2-colored-adam`.

``` {.bash}
  > mkdir ~/CGAL/builds/Triangulation_2-colored-adam
  > cd ~/CGAL/builds/Triangulation_2-colored-adam
  > cmake /path/to/my/cgal_repository
  > make
  > export CGAL_DIR=~/CGAL/builds/Triangulation_2-colored-adam
```

In such a setup she probably keeps a build of CGAL for each branch:

``` {.bash}
  > ls ~/CGAL/builds
  Convex_hull_4-jenny
  Convex_hull_2-make_it_faster-jenny
```

Jenny only has to set `$CGAL_DIR` to the build she wants to use before
building a demo

``` {.bash}
   > export CGAL_DIR=~/CGAL/builds/Convex_hull_2-make_it_faster-jenny
```

Building Adam's demo
====================

Now with properly set `$CGAL_DIR` it's time for Jenny to build a demo
contained from another branch in the local repository:

``` {.bash}
  > cd /path/to/my/cgal_repository/Triangulation_2/demo/Triangulation_2/
  > git branch
  master
  * Triangulation_2-coloring-adam
  ...
  > ls
  CMakeLists.txt colored_t2.cpp
```

It might be required to generate the file `CMakeLists.txt`. This can be
done with the script [`cgal_create_cmake_script`](Scripts/scripts/cgal_create_cmake_script)
found in `/path/to/my/cgal/Scripts/scripts/`. Then call:

``` {.bash}
  > cmake .
  > make
  > ./colored_t2
```

'''The important fact is that all headers are taken from the workingcopy
that have been used while building CGAL in the current `$CGAL_DIR`.

Remarks on this:

-   A build may be outdated after an include/dir has been deleted,
    switched or even updated. This might lead to compile problems (link
    with outdated version). Thus, it is recommended to build CGAL after
    each update, switch, merge of a branch (in particular if directories
    have been added/deleted, or cpp files have been added, deleted or
    altered).
-   There is currently no warning that a build does not match the
    branch. However, a little experience helps to be aware of this
    problem.
-   It might be possible to add such a warning (work in progress).

Finally, she's happy about Adam's cute pictures.

No installation required
========================

When you install CGAL (using `make install`), you copy all generated
library-objects and header files from the branch to a predefined
directory, `CMAKE_INSTALLATION_PREFIX`. Then, if you want to use the
installed version, you need to set
`$CGAL_DIR=$CMAKE_INSTALLATION_PREFIX` accordingly, thereby detaching
the installed files from the branch and loosing the connection to the
source-code management system. For most purposes there is no need to
install CGAL.
