# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( Surface_mesh_simplification_Examples )


cmake_minimum_required(VERSION 3.1)



# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# Boost and its components
find_package( Boost )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()


find_package( OpenMesh QUIET )

if ( OpenMesh_FOUND )
include( UseOpenMesh )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()

# include for local directory

# include for local package

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program( "edge_collapse_constrain_sharp_edges.cpp" )
create_single_source_cgal_program( "edge_collapse_constrained_border_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_enriched_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_polyhedron.cpp" )
create_single_source_cgal_program( "edge_collapse_surface_mesh.cpp" )
create_single_source_cgal_program( "edge_collapse_linear_cell_complex.cpp" )
create_single_source_cgal_program( "edge_collapse_constrained_border_surface_mesh.cpp" )
create_single_source_cgal_program( "edge_collapse_all_short_edges.cpp" )
create_single_source_cgal_program( "edge_collapse_bounded_normal_change.cpp" )
create_single_source_cgal_program( "edge_collapse_visitor_surface_mesh.cpp" )

if(OpenMesh_FOUND)
  create_single_source_cgal_program( "edge_collapse_OpenMesh.cpp" )

  target_link_libraries( edge_collapse_OpenMesh PRIVATE ${OPENMESH_LIBRARIES} )
endif()
