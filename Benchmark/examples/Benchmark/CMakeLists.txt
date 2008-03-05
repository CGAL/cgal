# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Benchmark_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( bbox.cpp )
  create_single_source_cgal_program( check_syntax.cpp )
  create_single_source_cgal_program( leftturn.cpp )
  create_single_source_cgal_program( polynomial_visitor.cpp )
  create_single_source_cgal_program( simple.cpp )

endif()
