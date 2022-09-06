# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

project(BGL_graphcut_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

create_single_source_cgal_program("alpha_expansion_example.cpp")
create_single_source_cgal_program("face_selection_borders_regularization_example.cpp")
