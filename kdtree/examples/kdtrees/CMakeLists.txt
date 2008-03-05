# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (kdtrees_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( example1.cpp )
  create_single_source_cgal_program( example2.cpp )
  create_single_source_cgal_program( example3.cpp )
  create_single_source_cgal_program( example4.cpp )

endif()
