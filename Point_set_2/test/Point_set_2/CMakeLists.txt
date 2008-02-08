# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( nearest_nb1.cpp )
  create_single_source_cgal_program( nearest_nb_fcn.cpp )
  create_single_source_cgal_program( range_search.cpp )
  create_single_source_cgal_program( range_search_fcn.cpp )
  create_single_source_cgal_program( rs_check_empty.cpp )

endif()
