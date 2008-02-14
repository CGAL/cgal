# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Algebraic_foundations_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( $ENV{CGAL_DIR} )
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "$ENV{CGAL_DIR}/cmake/modules"  )
endif()

find_package(CGAL)

if ( CGAL_FOUND )

  include(UseCGAL)
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( algebraic_structure_dispatch.cpp )
  create_single_source_cgal_program( fraction_traits.cpp )
  create_single_source_cgal_program( implicit_interoperable_dispatch.cpp )
  create_single_source_cgal_program( integralize.cpp )
  create_single_source_cgal_program( interoperable.cpp )

endif()
