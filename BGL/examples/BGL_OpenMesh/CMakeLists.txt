# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.29)
project(BGL_OpenMesh_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(OpenMesh)
if(OpenMesh_FOUND)
  include(CGAL_OpenMesh_support)
  create_single_source_cgal_program("TriMesh.cpp")
  target_link_libraries(TriMesh PRIVATE CGAL::OpenMesh_support)
else()
  message("NOTICE: This project requires OpenMesh and will not be compiled.")
endif()
