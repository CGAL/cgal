# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( approx_convex_partition_2_test_C.cpp )
  create_single_source_cgal_program( approx_convex_partition_2_test_H.cpp )
  create_single_source_cgal_program( greene_approx_convex_partition_2_test_C.cpp )
  create_single_source_cgal_program( greene_approx_convex_partition_2_test_H.cpp )
  create_single_source_cgal_program( is_y_monotone_2_test.cpp )
  create_single_source_cgal_program( optimal_convex_partition_2_test_C.cpp )
  create_single_source_cgal_program( optimal_convex_partition_2_test_H.cpp )
  create_single_source_cgal_program( partition_2_defaults_test_C.cpp )
  create_single_source_cgal_program( partition_2_defaults_test_H.cpp )
  create_single_source_cgal_program( y_monotone_partition_2_test_C.cpp )
  create_single_source_cgal_program( y_monotone_partition_2_test_H.cpp )

endif()
