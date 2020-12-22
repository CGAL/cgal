# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Triangulation_Examples)

if(CMAKE_COMPILER_IS_GNUCCX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
  message(
    STATUS
      "NOTICE: this directory requires a version of gcc >= 4.4, and will not be compiled."
  )
  return()
endif()

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)
  create_single_source_cgal_program("barycentric_subdivision.cpp")
  create_single_source_cgal_program("delaunay_triangulation.cpp")
  create_single_source_cgal_program("convex_hull.cpp")
  create_single_source_cgal_program("regular_triangulation.cpp")
  create_single_source_cgal_program("triangulation.cpp")
  create_single_source_cgal_program("triangulation_data_structure_dynamic.cpp")
  create_single_source_cgal_program("triangulation_data_structure_static.cpp")

  foreach(
    target
    barycentric_subdivision
    delaunay_triangulation
    convex_hull
    regular_triangulation
    triangulation
    triangulation_data_structure_dynamic
    triangulation_data_structure_static)
    target_link_libraries(${target} PUBLIC CGAL::Eigen_support)
  endforeach()
else()
  message(
    STATUS
      "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled."
  )
endif()
