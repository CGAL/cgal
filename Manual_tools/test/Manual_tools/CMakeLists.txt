# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Manual_tools_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)


endif()
