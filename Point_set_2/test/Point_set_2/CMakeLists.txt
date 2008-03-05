# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Point_set_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( nearest_nb1.cpp )
  create_single_source_cgal_program( nearest_nb_fcn.cpp )
  create_single_source_cgal_program( range_search.cpp )
  create_single_source_cgal_program( range_search_fcn.cpp )
  create_single_source_cgal_program( rs_check_empty.cpp )

endif()
