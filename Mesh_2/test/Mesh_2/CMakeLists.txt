# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Mesh_2_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()

if(cxx_std_20 IN_LIST CMAKE_CXX_COMPILE_FEATURES)
  target_compile_features(test_mesh_terrain PRIVATE cxx_std_20)
endif()
