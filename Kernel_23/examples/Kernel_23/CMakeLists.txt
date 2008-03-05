# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Kernel_23_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( MyKernel.cpp )
  create_single_source_cgal_program( test_delaunay_2.cpp )
  create_single_source_cgal_program( test_polyhedron_prog_incr_builder.cpp )
  create_single_source_cgal_program( test_polyhedron_prog_incr_builder_own_archetype.cpp )

endif()
