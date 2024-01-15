# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Alpha_wrap_3_Examples)

find_package(CGAL REQUIRED)

# create a target per cppfile
create_single_source_cgal_program("triangle_mesh_wrap.cpp")
create_single_source_cgal_program("triangle_soup_wrap.cpp")
create_single_source_cgal_program("point_set_wrap.cpp")
create_single_source_cgal_program("wrap_from_cavity.cpp")
create_single_source_cgal_program("mixed_inputs_wrap.cpp")
create_single_source_cgal_program("volumetric_wrap.cpp")
create_single_source_cgal_program("successive_wraps.cpp")
create_single_source_cgal_program("pause_and_resume_wrapping.cpp")
