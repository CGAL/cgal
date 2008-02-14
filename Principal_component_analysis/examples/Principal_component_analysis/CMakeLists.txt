# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Principal_component_analysis_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( barycenter.cpp )
  create_single_source_cgal_program( bounding_box.cpp )
  create_single_source_cgal_program( centroid.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_circles_2.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_cuboids_3.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_points_2.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_points_3.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_rectangles_2.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_segments_2.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_segments_3.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_spheres_3.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_tetrahedra_3.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_triangles_2.cpp )
  create_single_source_cgal_program( linear_least_squares_fitting_triangles_3.cpp )

endif()
