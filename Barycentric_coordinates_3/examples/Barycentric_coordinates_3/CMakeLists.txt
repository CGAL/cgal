# Created by the script cgal_create_cmake_script.
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Barycentric_coordinates_3_Examples)

find_package(CGAL REQUIRED)

include(CGAL_CreateSingleSourceCGALProgram)

create_single_source_cgal_program("tetrahedron_coordinates.cpp")
create_single_source_cgal_program("shape_deformation_3.cpp")
create_single_source_cgal_program("wachspress_coordinates_3.cpp")
create_single_source_cgal_program("discrete_harmonic_coordinates_3.cpp")
create_single_source_cgal_program("mean_value_coordinates_3.cpp")
