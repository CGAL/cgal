 S t r a i g h t S k e l
=========================

StraightSkel is an implementation of Straight Skeleton in 2- and 3-dimensional space.
It is used to compute offsets of closed polygons and polyhedrons.

The straight skeleton is always interior to the shape, and the offsetting is always inwards.

 Requirements
--------------

* C++ Compiler
  http://gcc.gnu.org/
* CMake build system
  http://www.cmake.org/
* Boost C++ Libraries
  http://www.boost.org/
* Computational Geometry Algorithms Library (CGAL)
  http://www.cgal.org/

Optional:
* SQLite
  http://www.sqlite.org/
* OpenGL Utility Toolkit (GLUT)
  http://freeglut.sourceforge.net/
* Doxygen
  http://www.stack.nl/~dimitri/doxygen/


 Building
----------

CMake is used as build system.
Therefore you may use "ccmake" or "cmake-gui" to configure the build process.

 Running
---------

The main executable is `StraightSkel`. Its usage is printed to standard output when the command
is invoked with no options.

Various advanced options can be found in the configuration file called "StraightSkel.ini".
This file is located in the same directory as the executable itself (build directory).

For outward offset mesh generation, a helper function is provided, which handles creating
an outside bounding box, inverting the mesh, etc. See sample code: src/offset_mesh.cpp

 Testing
---------

"test.sh" tests the code by generating skeleton & offset of many simple polyhedrons.
A summary is shown at the end.

Unit tests also document the usage of the code.
The build process needs to be configured to build the tests.
"make test" will sequentially run all "*TestRunner".
"*TestRunner --log_level=all" gives detailed output.

Debugging:
gdb --args ./StraightSkel 3d load anything.obj

Profiling:

valgrind --tool=callgrind ./StraightSkel 3d load anything.obj --no-window
kcachegrind callgrind.out.*

 Development
-------------

To avoid segmentation fault errors, standard C pointers are used very rarely.
A slightly modified version of the STL's Smart Pointer may be used instead.
(Modification: Print the stack trace in case of an error.)

Organization of Source Code:

src/
  algo/     Algorithms
  data/     Data Structures
  db/       Database
  kernel/   Geometric Kernel (double precision)
  ui/       User Interface
    gl/     OpenGL
    ps/     PostScript
  util/     Tools
  misc/     Miscellaneous (pre- and post-processing helpers for mesh offsetting)
test/       Unit Test Cases

 Database
----------

sqlite3 skeldata2d.db3 "SELECT * FROM Polygons"
sqlite3 skeldata3d.db3 "SELECT * FROM Polyhedrons"

 Visualization
---------------

Keybindings for the OpenGL window can be found and modified in the configuration file ("StraightSkel.ini",
in the build directory).

Creating animations:
The OpenGL window can be dumped to a bitmap file (.bmp) by pressing the associated button. (default: F10)
A series of bitmap files may be converted to an animated GIF using ImageMagick's convert:

```
convert  -delay 100  -loop 0  -colors 64  *.bmp  animated.gif
```

To create a movie you may use mencoder:

```
mencoder 'mf://*.bmp' -mf type=bmp:fps=10 -ovc xvid -xvidencopts bitrate=1000 -o output.avi
```

 History
---------

StraightSkel was started by Gernot Walzl in the years 2011, 2012, 2013 (MIT).

StraightSkel was further developed by GeometryFactory since 2024 (GPL-3.0-or-later OR LicenseRef-Commercial).


 Disclaimer (MIT code)
----------------------

THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
