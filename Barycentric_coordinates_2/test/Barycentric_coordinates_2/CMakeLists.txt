# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)

project(Barycentric_coordinates_2_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

create_single_source_cgal_program("test_almost_degenerate_segment.cpp")
create_single_source_cgal_program("test_segment_coordinates.cpp")
create_single_source_cgal_program("test_segment_coordinates_with_offset.cpp")

create_single_source_cgal_program("test_almost_degenerate_triangle.cpp")
create_single_source_cgal_program("test_triangle_coordinates.cpp")

create_single_source_cgal_program("test_boundary_coordinates_at_vertices.cpp")
create_single_source_cgal_program("test_boundary_coordinates_on_edges.cpp")
create_single_source_cgal_program("test_bc_cgal_polygons.cpp")

create_single_source_cgal_program("test_wp_dh_unit_square.cpp")
create_single_source_cgal_program("test_wp_almost_degenerate_polygon.cpp")
create_single_source_cgal_program("test_wp_const_linear_precision.cpp")
create_single_source_cgal_program("test_wp_triangle.cpp")
create_single_source_cgal_program("test_wp_weights.cpp")

create_single_source_cgal_program("test_mv_special_points.cpp")
create_single_source_cgal_program("test_mv_weakly_convex_polygon.cpp")
create_single_source_cgal_program("test_mv_const_linear_precision.cpp")
create_single_source_cgal_program("test_mv_triangle.cpp")
create_single_source_cgal_program("test_mv_weights.cpp")

create_single_source_cgal_program("test_dh_almost_degenerate_polygon.cpp")
create_single_source_cgal_program("test_dh_const_linear_precision.cpp")
create_single_source_cgal_program("test_dh_triangle.cpp")
create_single_source_cgal_program("test_dh_weights.cpp")

# this code is deprecated:
create_single_source_cgal_program("test_sc_deprecated_api.cpp")
create_single_source_cgal_program("test_tc_deprecated_api.cpp")
create_single_source_cgal_program("test_wp_deprecated_api.cpp")
create_single_source_cgal_program("test_mv_deprecated_api.cpp")
create_single_source_cgal_program("test_dh_deprecated_api.cpp")

find_package(Eigen3 3.1.0 QUIET) # (3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("test_hm_unit_square.cpp")
  target_link_libraries(test_hm_unit_square PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("test_hm_const_linear_precision.cpp")
  target_link_libraries(test_hm_const_linear_precision PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("test_hm_triangle.cpp")
  target_link_libraries(test_hm_triangle PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("test_bc_projection_traits.cpp")
  target_link_libraries(test_bc_projection_traits PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("test_bc_all_coordinates.cpp")
  target_link_libraries(test_bc_all_coordinates PUBLIC CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: Several tests require the Eigen library, and will not be compiled.")
endif()
