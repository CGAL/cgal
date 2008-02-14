# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Interpolation_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( interpolation_2.cpp )
  create_single_source_cgal_program( linear_interpolation_2.cpp )
  create_single_source_cgal_program( nn_coordinates_2.cpp )
  create_single_source_cgal_program( nn_coordinates_3.cpp )
  create_single_source_cgal_program( rn_coordinates_2.cpp )
  create_single_source_cgal_program( sibson_interpolation_2.cpp )
  create_single_source_cgal_program( surface_neighbor_coordinates_3.cpp )

endif()
