# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(Spatial_searching_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("circular_query.cpp")
create_single_source_cgal_program("distance_browsing.cpp")
create_single_source_cgal_program("iso_rectangle_2_query.cpp")
create_single_source_cgal_program("iso_rectangle_2_query_projection.cpp")
create_single_source_cgal_program("nearest_neighbor_searching.cpp")
create_single_source_cgal_program("searching_with_circular_query.cpp")
create_single_source_cgal_program("searching_with_point_with_info.cpp")
create_single_source_cgal_program("searching_with_point_with_info_inplace.cpp")
create_single_source_cgal_program("searching_with_point_with_info_pmap.cpp")
create_single_source_cgal_program("searching_surface_mesh_vertices.cpp")
create_single_source_cgal_program("searching_polyhedron_vertices.cpp")
create_single_source_cgal_program("searching_triangulation_vertices.cpp")
create_single_source_cgal_program("searching_polyhedron_vertices_with_fuzzy_sphere.cpp")
create_single_source_cgal_program("user_defined_point_and_distance.cpp")
create_single_source_cgal_program("using_fair_splitting_rule.cpp")
create_single_source_cgal_program("weighted_Minkowski_distance.cpp")
create_single_source_cgal_program("searching_sphere_orthogonally.cpp")
create_single_source_cgal_program("splitter_worst_cases.cpp")

find_package(Eigen3 3.1.91 QUIET) #(requires 3.1.91 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("fuzzy_range_query.cpp")
  target_link_libraries(fuzzy_range_query PUBLIC CGAL::Eigen3_support)

  create_single_source_cgal_program("general_neighbor_searching.cpp")
  target_link_libraries(general_neighbor_searching PUBLIC CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: The examples 'fuzzy_range_query' and 'general_neighbor_searching'")
  message(STATUS "will not be compiled as they use CGAL::Epick_d, which requires the Eigen library.")
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  message(STATUS "Found TBB")

  create_single_source_cgal_program("parallel_kdtree.cpp")
  target_link_libraries(parallel_kdtree PUBLIC CGAL::TBB_support)
else()
  message(STATUS "NOTICE: The example 'parallel_kdtree' requires TBB, and will not be compiled")
endif()
