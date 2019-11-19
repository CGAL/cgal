# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)
project( BGL_Tests )





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
  add_definitions( -DCGAL_USE_OPENMESH )
else()
  message(STATUS "Examples that use OpenMesh will not be compiled.")
endif()



# include for local package

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

if(OpenMesh_FOUND)
  create_single_source_cgal_program( "graph_concept_OpenMesh.cpp" )
  target_link_libraries( graph_concept_OpenMesh PRIVATE ${OPENMESH_LIBRARIES} )
endif()

create_single_source_cgal_program( "next.cpp" )

create_single_source_cgal_program( "test_circulator.cpp" )

create_single_source_cgal_program( "test_Gwdwg.cpp" )

create_single_source_cgal_program( "test_bgl_dual.cpp" )

create_single_source_cgal_program( "graph_concept_Polyhedron_3.cpp" )

create_single_source_cgal_program( "graph_concept_Dual.cpp" )

create_single_source_cgal_program( "graph_concept_Triangulation_2.cpp" )

create_single_source_cgal_program( "graph_concept_Surface_mesh.cpp" )

create_single_source_cgal_program( "graph_concept_Seam_mesh_Surface_mesh.cpp" )

create_single_source_cgal_program( "graph_concept_Gwdwg_Surface_mesh.cpp" )

create_single_source_cgal_program( "graph_concept_Linear_cell_complex.cpp" )

create_single_source_cgal_program( "graph_concept_Arrangement_2.cpp" )

create_single_source_cgal_program( "test_clear.cpp" )

create_single_source_cgal_program( "test_helpers.cpp" )

create_single_source_cgal_program( "test_Has_member_clear.cpp" )

create_single_source_cgal_program( "test_Has_member_id.cpp" )

create_single_source_cgal_program( "test_cgal_bgl_named_params.cpp" )

create_single_source_cgal_program( "test_bgl_read_write.cpp" )

create_single_source_cgal_program( "graph_concept_Face_filtered_graph.cpp" )

create_single_source_cgal_program( "test_Manifold_face_removal.cpp")

create_single_source_cgal_program( "test_Face_filtered_graph.cpp" )

create_single_source_cgal_program( "test_Euler_operations.cpp" )

create_single_source_cgal_program( "test_Collapse_edge.cpp" )

create_single_source_cgal_program( "test_graph_traits.cpp" )

create_single_source_cgal_program( "test_Properties.cpp" )

create_single_source_cgal_program( "test_wrl.cpp" )

if(OpenMesh_FOUND)
  target_link_libraries( test_clear PRIVATE ${OPENMESH_LIBRARIES})
  target_link_libraries( test_Euler_operations PRIVATE ${OPENMESH_LIBRARIES})
  target_link_libraries( test_Collapse_edge PRIVATE ${OPENMESH_LIBRARIES})
  target_link_libraries( test_Face_filtered_graph PRIVATE ${OPENMESH_LIBRARIES})
  target_link_libraries( test_graph_traits PRIVATE ${OPENMESH_LIBRARIES} )
  target_link_libraries( test_Properties PRIVATE ${OPENMESH_LIBRARIES})
  target_link_libraries( test_bgl_read_write PRIVATE ${OPENMESH_LIBRARIES})
endif()
