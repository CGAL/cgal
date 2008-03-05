# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Spatial_searching_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( Building_kd_tree_with_ddim_points.cpp )
  create_single_source_cgal_program( Building_kd_tree_with_own_pointtype.cpp )
  create_single_source_cgal_program( Circular_query.cpp )
  create_single_source_cgal_program( Iso_rectangle_2_query.cpp )
  create_single_source_cgal_program( K_neighbor_search_manhattan_distance_isobox_point.cpp )
  create_single_source_cgal_program( K_neighbor_search_with_circle.cpp )
  create_single_source_cgal_program( Orthogonal_incremental_neighbor_search.cpp )
  create_single_source_cgal_program( Orthogonal_k_neighbor_search.cpp )
  create_single_source_cgal_program( Range_searching.cpp )
  create_single_source_cgal_program( Splitters.cpp )

endif()
