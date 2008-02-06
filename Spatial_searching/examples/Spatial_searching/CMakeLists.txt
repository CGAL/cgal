# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( circular_query.cpp )
  create_single_source_cgal_program( distance_browsing.cpp )
  create_single_source_cgal_program( fuzzy_range_query.cpp )
  create_single_source_cgal_program( general_neighbor_searching.cpp )
  create_single_source_cgal_program( iso_rectangle_2_query.cpp )
  create_single_source_cgal_program( nearest_neighbor_searching.cpp )
  create_single_source_cgal_program( searching_with_circular_query.cpp )
  create_single_source_cgal_program( user_defined_point_and_distance.cpp )
  create_single_source_cgal_program( using_fair_splitting_rule.cpp )
  create_single_source_cgal_program( weighted_Minkowski_distance.cpp )

endif()
