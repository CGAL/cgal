# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Circular_kernel_3_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( spherical_kernel_3.cpp )

endif()
