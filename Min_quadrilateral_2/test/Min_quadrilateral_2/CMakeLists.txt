# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGALROOT}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( minimum_enclosing_quadrilateral_2_test_C.cpp )
  create_single_source_cgal_program( minimum_enclosing_quadrilateral_2_test_H.cpp )
  create_single_source_cgal_program( minimum_enclosing_quadrilateral_2_test_traits.cpp )

endif()
