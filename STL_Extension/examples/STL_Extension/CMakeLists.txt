# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (STL_Extension_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( in_place_list_prog.cpp )
  create_single_source_cgal_program( min_element_if_example_noheader.cpp )
  create_single_source_cgal_program( min_max_element_example_noheader.cpp )

endif()
