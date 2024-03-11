 S t r a i g h t S k e l
=========================

StraightSkel is an implementation of the Straight Skeleton in 2- and
3-dimensional space. It is used to animate the computation of offsets
of polygons and polyhedrons.
StraightSkel was written by Gernot Walzl in the years 2011, 2012, 2013.


 Requirements
--------------

* C++ Compiler  (of course)
  http://gcc.gnu.org/
* Boost C++ Libraries
  http://www.boost.org/
* CMake build system
  http://www.cmake.org/
* SQLite
  http://www.sqlite.org/
* OpenGL Utility Toolkit (GLUT)  (often part of the graphics card driver)
  http://freeglut.sourceforge.net/

Optional:
* Computational Geometry Algorithms Library (CGAL)
  http://www.cgal.org/
* Doxygen
  http://www.stack.nl/~dimitri/doxygen/


 Building
----------

"build.sh" does an "out-of-source build".
It creates the binaries in a directory called "build".
A database of some given polygons is created also.

CMake is used as build system.
Therefore you may use "ccmake" or "cmake-gui" to configure the build process.

This is useful to switch to the exact geometric kernel of CGAL.
The included kernel uses double precision only.


 Running
---------

Key bindings and various options may be found in the configuration file
called "StraightSkel.ini". This file is located in the same directory as the
executable itself.

The usage of StraightSkel is printed to standard output when the command
is invoked with no options.

Common problems:
* Plane equations are not exact:
  Many commonly used file formats (like Wavefront's .obj) have a
  number representation with limitations. It is not always possible
  that more than 3 points are exactly on the same plane.
* Degenerated input:
  More than 1 event occurs exactly on the same offset.

Creating animations:
The OpenGL window can be dumped to a bitmap file (.bmp) by pressing the
associated button. (default: F10)
A series of bitmap files may be converted to an animated GIF using
ImageMagick's convert:

convert  -delay 100  -loop 0  -colors 64  *.bmp  animated.gif

To create a movie you may use mencoder:

mencoder 'mf://*.bmp' -mf type=bmp:fps=10 -ovc xvid -xvidencopts bitrate=1000 -o output.avi


 Documentation
---------------

The build process can be configured to use doxygen to create documentation
from the source code.
Once the makefile is configured, "make doc" will invoke doxygen.

Additional documentation may be found inside the "doc" folder.


 Development
-------------

To avoid segmentation fault errors, standard C pointers are used very rarely.
A slightly modified version of Boost's Smart Pointer is used instead.
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
  util/
test/       Unit Test Cases

Testing:
Testing is a very important part to improve the quality of the code.
"test.sh" tests the executable by invoking the command on many polyhedrons.
A summary is shown at the end.

Unit Test Cases:
Unit tests also document the usage of the code.
The build process needs to be configured to build the tests.
"make test" will sequentially run all "*TestRunner".
"*TestRunner --log_level=all" gives detailed output.

Debugging:
gdb is useful to detect errors.

gdb ./StraightSkel
run 3d load anything.obj

Profiling:

valgrind --tool=callgrind ./StraightSkel 3d load anything.obj --no-window
kcachegrind callgrind.out.*


 Database
----------

sqlite3 skeldata2d.db3 "SELECT * FROM Polygons"
sqlite3 skeldata3d.db3 "SELECT * FROM Polyhedrons"


 Disclaimer
------------

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
