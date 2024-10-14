# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Surface_mesh_simplification_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("edge_collapse_envelope.cpp")
create_single_source_cgal_program("edge_collapse_constrain_sharp_edges.cpp")
create_single_source_cgal_program("edge_collapse_constrained_border_polyhedron.cpp")
create_single_source_cgal_program("edge_collapse_enriched_polyhedron.cpp")
create_single_source_cgal_program("edge_collapse_polyhedron.cpp")
create_single_source_cgal_program("edge_collapse_surface_mesh.cpp")
create_single_source_cgal_program("edge_collapse_linear_cell_complex.cpp")
create_single_source_cgal_program("edge_collapse_constrained_border_surface_mesh.cpp")
create_single_source_cgal_program("edge_collapse_all_short_edges.cpp")
create_single_source_cgal_program("edge_collapse_bounded_normal_change.cpp")
create_single_source_cgal_program("edge_collapse_visitor_surface_mesh.cpp")

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("edge_collapse_garland_heckbert.cpp")
  target_link_libraries(edge_collapse_garland_heckbert PUBLIC CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Garland-Heckbert polices require the Eigen library, which has not been found; related examples will not be compiled.")
endif()

find_package(OpenMesh QUIET)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)
  message(STATUS "Found OpenMesh")

  create_single_source_cgal_program("edge_collapse_OpenMesh.cpp")
  target_link_libraries(edge_collapse_OpenMesh PRIVATE CGAL::OpenMesh_support)
else()
  message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
endif()

find_package(METIS QUIET)
include(CGAL_METIS_support)

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support AND TARGET CGAL::METIS_support)
  message(STATUS "Found METIS & TBB")

  create_single_source_cgal_program("collapse_small_edges_in_parallel.cpp")
  target_link_libraries(collapse_small_edges_in_parallel PUBLIC CGAL::TBB_support CGAL::METIS_support)
else()
  message(STATUS "NOTICE: The example 'collapse_small_edges_in_parallel' uses the METIS and TBB libraries, and will not be compiled.")
endif()
