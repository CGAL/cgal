# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(PMP_Boolean_operations_Benchmarks)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("rotated_cubes_autorefinement.cpp")
create_single_source_cgal_program("coplanar_cubes_autorefinement.cpp")

create_single_source_cgal_program("Performance/performance_snap_polygon_soup.cpp")
create_single_source_cgal_program("Robustness/robustness_snap_polygon_soup.cpp")
create_single_source_cgal_program("Quality/quality_snap_polygon_soup.cpp")

create_single_source_cgal_program("polygon_mesh_slicer.cpp")
target_link_libraries(polygon_mesh_slicer PRIVATE CGAL::Eigen3_support)

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(rotated_cubes_autorefinement PRIVATE CGAL::TBB_support)
  target_link_libraries(coplanar_cubes_autorefinement PRIVATE CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Sequential code will be used.")
endif()
