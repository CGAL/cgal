# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)
project(Hash_map)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

create_single_source_cgal_program("hm.cpp")
create_single_source_cgal_program("foreach.cpp")
create_single_source_cgal_program("triangulation.cpp")
create_single_source_cgal_program("polyhedron.cpp")
create_single_source_cgal_program("surface_mesh.cpp")
create_single_source_cgal_program("arrangement.cpp")
