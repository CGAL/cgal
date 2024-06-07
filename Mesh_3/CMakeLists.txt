# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.1...3.23)
project(Mesh_3)

if(CGAL_ENABLE_TESTING)
  enable_testing()
endif()

option(CGAL_MESH_3_BENCHMARK "Add benchmark to project" OFF)
if(CGAL_MESH_3_BENCHMARK)
  add_subdirectory("benchmark/Mesh_3")
endif()

add_subdirectory("examples/Mesh_3")
add_subdirectory("test/Mesh_3")
