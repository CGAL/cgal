# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Min_quadrilateral_2_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( minimum_enclosing_quadrilateral_2_test_C.cpp )
  create_single_source_cgal_program( minimum_enclosing_quadrilateral_2_test_H.cpp )
  create_single_source_cgal_program( minimum_enclosing_quadrilateral_2_test_traits.cpp )

endif()
