# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( random_convex_set.cpp )
  create_single_source_cgal_program( random_degenerate_point_set.cpp )
  create_single_source_cgal_program( random_grid.cpp )
  create_single_source_cgal_program( random_polygon2.cpp )
  create_single_source_cgal_program( random_polygon.cpp )
  create_single_source_cgal_program( random_segments1.cpp )
  create_single_source_cgal_program( random_segments2.cpp )

endif()
