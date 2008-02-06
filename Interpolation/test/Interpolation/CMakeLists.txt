# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( test_interpolation_functions_2.cpp )
  create_single_source_cgal_program( test_natural_neighbors_2.cpp )
  create_single_source_cgal_program( test_regular_neighbors_2.cpp )
  create_single_source_cgal_program( test_surface_neighbor_3.cpp )

endif()
