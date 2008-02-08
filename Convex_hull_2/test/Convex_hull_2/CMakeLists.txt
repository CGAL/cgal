# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  include_directories (BEFORE include)

  create_single_source_cgal_program( ch_test_CH.cpp )
  create_single_source_cgal_program( ch_test_SC.cpp )
  create_single_source_cgal_program( ch_test_SH.cpp )
  create_single_source_cgal_program( ch_test_SS.cpp )
  create_single_source_cgal_program( ka_test_convex_hull_2.cpp )

endif()
