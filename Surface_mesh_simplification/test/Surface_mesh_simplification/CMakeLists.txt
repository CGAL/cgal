# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Surface_mesh_simplification_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( test_edge_collapse_Polyhedron_3.cpp )
  create_single_source_cgal_program( test_LT_strategy.cpp )
  create_single_source_cgal_program( test_MP_strategy.cpp )

endif()
