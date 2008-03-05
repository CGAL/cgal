# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Algebraic_foundations_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( algebraic_structure_dispatch.cpp )
  create_single_source_cgal_program( fraction_traits.cpp )
  create_single_source_cgal_program( implicit_interoperable_dispatch.cpp )
  create_single_source_cgal_program( integralize.cpp )
  create_single_source_cgal_program( interoperable.cpp )

endif()
