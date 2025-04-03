# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Barycentric_coordinates_3_Tests)

find_package(CGAL REQUIRED)

include(CGAL_CreateSingleSourceCGALProgram)

create_single_source_cgal_program("test_tetrahedron_coordinates.cpp")
create_single_source_cgal_program("test_wp_dh_mv_tetrahedron.cpp")
create_single_source_cgal_program("test_wp_weights_3.cpp")
create_single_source_cgal_program("test_dh_weights_3.cpp")
create_single_source_cgal_program("test_mv_weights_3.cpp")
create_single_source_cgal_program("test_edge_cases.cpp")
create_single_source_cgal_program("test_containers.cpp")
create_single_source_cgal_program("test_boundary_coordinates.cpp")
