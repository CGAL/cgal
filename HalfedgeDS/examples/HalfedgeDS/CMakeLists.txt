# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (HalfedgeDS_Example)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})
  include(CreateSingleSourceCGALProgram)

  create_single_source_cgal_program( hds_prog_color.cpp )
  create_single_source_cgal_program( hds_prog_compact2.cpp )
  create_single_source_cgal_program( hds_prog_compact.cpp )
  create_single_source_cgal_program( hds_prog_default.cpp )
  create_single_source_cgal_program( hds_prog_default_with_ccb.cpp )
  create_single_source_cgal_program( hds_prog_edge_iterator.cpp )
  create_single_source_cgal_program( hds_prog_graph2.cpp )
  create_single_source_cgal_program( hds_prog_graph.cpp )
  create_single_source_cgal_program( hds_prog_graph_with_ccb.cpp )
  create_single_source_cgal_program( hds_prog_halfedge_iterator.cpp )
  create_single_source_cgal_program( hds_prog_vector.cpp )

endif()
