# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Box_intersection_d_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( automated_test.cpp )
  create_single_source_cgal_program( benchmark.cpp )
  create_single_source_cgal_program( box_grid.cpp )
  create_single_source_cgal_program( random_set_test.cpp )

endif()
