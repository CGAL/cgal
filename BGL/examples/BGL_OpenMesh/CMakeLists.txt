# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)
project(BGL_OpenMesh_Examples)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(OpenMesh)
if(OpenMesh_FOUND)
  include(UseOpenMesh)
  create_single_source_cgal_program("TriMesh.cpp")
  target_link_libraries(TriMesh PRIVATE ${OPENMESH_LIBRARIES})
else()
  message("NOTICE: This project requires OpenMesh and will not be compiled.")
endif()
