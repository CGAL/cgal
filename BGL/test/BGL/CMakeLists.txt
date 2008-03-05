# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (BGL_Test)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( halfedge_graph_traits_Polyhedron_3.cpp )
  create_single_source_cgal_program( properties_Polyhedron_3.cpp )

endif()
