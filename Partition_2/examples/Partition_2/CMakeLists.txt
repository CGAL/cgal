# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Partition_2_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( approx_convex_partition_2.cpp )
  create_single_source_cgal_program( greene_approx_convex_partition_2.cpp )
  create_single_source_cgal_program( optimal_convex_partition_2.cpp )
  create_single_source_cgal_program( y_monotone_partition_2.cpp )

endif()
