# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(BGL_Tests)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("test_split.cpp")
create_single_source_cgal_program("next.cpp")
create_single_source_cgal_program("test_circulator.cpp")
create_single_source_cgal_program("test_Gwdwg.cpp")
create_single_source_cgal_program("test_bgl_dual.cpp")
create_single_source_cgal_program("graph_concept_Polyhedron_3.cpp")
create_single_source_cgal_program("graph_concept_Dual.cpp")
create_single_source_cgal_program("graph_concept_Triangulation_2.cpp")
create_single_source_cgal_program("graph_concept_Surface_mesh.cpp")
create_single_source_cgal_program("graph_concept_Seam_mesh_Surface_mesh.cpp")
create_single_source_cgal_program("graph_concept_Gwdwg_Surface_mesh.cpp")
create_single_source_cgal_program("graph_concept_Linear_cell_complex.cpp")
create_single_source_cgal_program("graph_concept_Arrangement_2.cpp")
create_single_source_cgal_program("graph_concept_Derived.cpp" )
create_single_source_cgal_program("test_clear.cpp" )
create_single_source_cgal_program("test_helpers.cpp")
create_single_source_cgal_program("test_Has_member_clear.cpp")
create_single_source_cgal_program("test_Has_member_id.cpp")
create_single_source_cgal_program("test_bgl_read_write.cpp")
create_single_source_cgal_program("graph_concept_Face_filtered_graph.cpp")
create_single_source_cgal_program("test_Manifold_face_removal.cpp")
create_single_source_cgal_program("test_Regularize_face_selection_borders.cpp")
create_single_source_cgal_program("test_Face_filtered_graph.cpp")
create_single_source_cgal_program("test_Euler_operations.cpp")
create_single_source_cgal_program("test_test_face.cpp" )
create_single_source_cgal_program("test_Collapse_edge.cpp" )
create_single_source_cgal_program("test_Collapse_edge_with_constraints.cpp" )
create_single_source_cgal_program("test_graph_traits.cpp")
create_single_source_cgal_program("test_Properties.cpp")
create_single_source_cgal_program("bench_read_from_stream_vs_add_face_and_add_faces.cpp")
create_single_source_cgal_program("graph_traits_inheritance.cpp" )
create_single_source_cgal_program("test_deprecated_io.cpp")

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)

  target_link_libraries(test_clear PRIVATE CGAL::OpenMesh_support)
  target_link_libraries(test_Euler_operations PRIVATE CGAL::OpenMesh_support)
  target_link_libraries(test_Collapse_edge PRIVATE CGAL::OpenMesh_support)
  target_link_libraries(test_Face_filtered_graph PRIVATE CGAL::OpenMesh_support)
  target_link_libraries(test_graph_traits PRIVATE CGAL::OpenMesh_support )
  target_link_libraries(test_Properties PRIVATE CGAL::OpenMesh_support)
  target_link_libraries(test_bgl_read_write PRIVATE CGAL::OpenMesh_support)

  create_single_source_cgal_program("graph_concept_OpenMesh.cpp")
  target_link_libraries(graph_concept_OpenMesh PRIVATE CGAL::OpenMesh_support)

  create_single_source_cgal_program("test_OpenMesh.cpp")
  target_link_libraries(test_OpenMesh PRIVATE CGAL::OpenMesh_support)


else()
  message(STATUS "NOTICE: Tests that use OpenMesh will not be compiled.")
endif()

find_package(VTK 9.0 QUIET COMPONENTS CommonCore IOCore IOLegacy IOXML FiltersCore FiltersSources)
if (VTK_FOUND AND VTK_LIBRARIES)
  message(STATUS "VTK ${VTK_VERSION} found ${VTK_LIBRARIES}")
  target_link_libraries(test_bgl_read_write PRIVATE ${VTK_LIBRARIES})
  target_compile_definitions(test_bgl_read_write PRIVATE -DCGAL_USE_VTK -DNOMINMAX)
  target_link_libraries(test_deprecated_io PRIVATE ${VTK_LIBRARIES})
  target_compile_definitions(test_deprecated_io PRIVATE -DCGAL_USE_VTK -DNOMINMAX)
else()
  message(STATUS "Tests that use VTK will not be compiled.")
endif() #VTK_FOUND

find_path(3MF_INCLUDE_DIR
  NAMES Model/COM/NMR_DLLInterfaces.h
  DOC "Path to lib3MF headers"
  )
find_library(3MF_LIBRARIES NAMES 3MF DOC "Path to the lib3MF library")
if(3MF_LIBRARIES AND 3MF_INCLUDE_DIR  AND EXISTS "${3MF_INCLUDE_DIR}/Model/COM/NMR_DLLInterfaces.h")
  include_directories(${3MF_INCLUDE_DIR})
  create_single_source_cgal_program( test_3mf_to_sm.cpp )
  target_link_libraries(test_3mf_to_sm PRIVATE ${3MF_LIBRARIES})
  target_compile_definitions(test_3mf_to_sm PRIVATE -DCGAL_LINKED_WITH_3MF)
else()
  message(STATUS "NOTICE: The test 'test_3mf_to_sm' requires the lib3MF library, and will not be compiled.")
endif()
