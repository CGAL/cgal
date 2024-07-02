# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Surface_mesh_simplification_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("edge_collapse_topology.cpp")
create_single_source_cgal_program("test_edge_collapse_bounded_distance.cpp")
create_single_source_cgal_program("test_edge_collapse_Envelope.cpp")
create_single_source_cgal_program("test_edge_collapse_Polyhedron_3.cpp")
create_single_source_cgal_program("test_edge_profile_link.cpp")
create_single_source_cgal_program("test_edge_deprecated_stop_predicates.cpp")
create_single_source_cgal_program("test_edge_collapse_stability.cpp")

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("edge_collapse_garland_heckbert_variations.cpp")
  target_link_libraries(edge_collapse_garland_heckbert_variations PUBLIC CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Garland-Heckbert polices require the Eigen library, which has not been found; related examples will not be compiled.")
endif()
