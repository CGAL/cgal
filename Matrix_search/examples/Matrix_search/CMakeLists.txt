# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Matrix_search_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( all_furthest_neighbors_2.cpp )
  create_single_source_cgal_program( extremal_polygon_2_area.cpp )
  create_single_source_cgal_program( extremal_polygon_2_perimeter.cpp )
  create_single_source_cgal_program( rectangular_p_center_2.cpp )
  create_single_source_cgal_program( sorted_matrix_search.cpp )

endif()
