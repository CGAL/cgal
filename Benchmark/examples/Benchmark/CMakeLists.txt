# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Benchmark_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( bbox.cpp )
  create_single_source_cgal_program( check_syntax.cpp )
  create_single_source_cgal_program( leftturn.cpp )
  create_single_source_cgal_program( polynomial_visitor.cpp )
  create_single_source_cgal_program( simple.cpp )

endif()
