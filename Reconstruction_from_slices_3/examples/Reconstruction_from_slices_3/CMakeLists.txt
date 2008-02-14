# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Reconstruction_from_slices_3_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( test_reconstruction_from_slices.cpp )
  create_single_source_cgal_program( test_triangulation_points_in_polylines.cpp )
  create_single_source_cgal_program( test_triangulation_points_in_slices.cpp )

endif()
