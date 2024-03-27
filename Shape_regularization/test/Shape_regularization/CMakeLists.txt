# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(Shape_regularization_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

# Find OSQP library and headers.
find_package(OSQP QUIET)
include(CGAL_OSQP_support)
if(TARGET CGAL::OSQP_support)
  message(STATUS "Found OSQP")

  set(osqp_targets
  test_segments
  test_overloads
  test_0_segments
  test_1_segment
  test_2_segments
  test_3_segments
  test_4_segments
  test_100_segments_angles
  test_100_segments_offsets
  test_cgal_solver)

  foreach(osqp_target ${osqp_targets})
    create_single_source_cgal_program("${osqp_target}.cpp")
    if(TARGET ${osqp_target})
      target_link_libraries(${osqp_target} PUBLIC CGAL::OSQP_support)
    endif()
  endforeach()
else()
  message(STATUS "NOTICE: OSQP was not found. OSQP tests won't be available.")
endif()

set(targets
test_neighbor_query
test_groups_parallel
test_groups_collinear
test_groups_orthogonal
test_unique_segments
test_plane_regularization
test_closed_contour_3_edges
test_closed_contour_4_edges
test_closed_contour_6_edges
test_closed_contour_10_edges
test_open_contour_1_edge
test_open_contour_2_edges
test_open_contour_3_edges
test_open_contour_5_edges
test_open_contour_9_edges
test_equal_contours
test_directions_longest
test_directions_multiple_1
test_directions_multiple_2
test_directions_multiple_3
test_directions_user_defined
test_contour_orientation)

foreach(target ${targets})
  create_single_source_cgal_program("${target}.cpp")
endforeach()
