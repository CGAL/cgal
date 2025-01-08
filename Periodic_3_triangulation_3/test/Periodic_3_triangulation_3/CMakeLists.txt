# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Periodic_3_triangulation_3_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

include_directories(BEFORE "include")

create_single_source_cgal_program("test_p3rt3_as_p3t3.cpp")
create_single_source_cgal_program("test_p3rt3_insert_remove_point_set.cpp")
create_single_source_cgal_program("test_p3rt3_versus_rt3.cpp")
create_single_source_cgal_program("test_periodic_3_alpha_shape_3.cpp")
create_single_source_cgal_program("test_periodic_3_delaunay_3.cpp")
create_single_source_cgal_program("test_periodic_3_delaunay_hierarchy_3.cpp")
create_single_source_cgal_program("test_periodic_3_offset_3.cpp")
create_single_source_cgal_program("test_periodic_3_regular_triangulation_3.cpp")
create_single_source_cgal_program(
  "test_periodic_3_regular_triangulation_traits_3.cpp")
create_single_source_cgal_program("test_periodic_3_triangulation_3.cpp")
create_single_source_cgal_program("test_periodic_3_triangulation_tds.cpp")
create_single_source_cgal_program(
  "test_periodic_3_triangulation_traits_C_3.cpp")
create_single_source_cgal_program(
  "test_periodic_3_triangulation_traits_H_3.cpp")
create_single_source_cgal_program(
  "test_periodic_3_triangulation_traits_SC_3.cpp")
create_single_source_cgal_program(
  "test_periodic_3_triangulation_traits_SH_3.cpp")
create_single_source_cgal_program("test_dummy_point_generation.cpp")
