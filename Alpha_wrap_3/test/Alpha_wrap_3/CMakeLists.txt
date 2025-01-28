# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Alpha_wrap_3_Tests)

find_package(CGAL REQUIRED)

# create a target per cppfile
create_single_source_cgal_program("test_alpha_wrap_3_mesh.cpp")
create_single_source_cgal_program("test_AW3_cavity_initializations.cpp")
create_single_source_cgal_program("test_AW3_manifoldness.cpp")
create_single_source_cgal_program("test_AW3_multiple_calls.cpp")
create_single_source_cgal_program("test_AW3_compilation.cpp")
