# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Set_movable_separability_2_Examples)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("top_edges_single_mold_trans_cast.cpp")
create_single_source_cgal_program("is_pullout_direction_single_mold_trans_cast.cpp")
create_single_source_cgal_program("pullout_directions_single_mold_trans_cast.cpp")
