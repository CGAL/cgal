# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Generator_Examples)

find_package(CGAL REQUIRED)

# Use Eigen
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen_support)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  if(NOT (${cppfile} STREQUAL "random_points_in_tetrahedral_mesh_3.cpp")
     OR NOT (${cppfile} STREQUAL "random_points_on_tetrahedral_mesh_3.cpp")
     OR TARGET CGAL::Eigen_support)
    create_single_source_cgal_program("${cppfile}")
    if(TARGET CGAL::Eigen_support)
      get_filename_component(target ${cppfile} NAME_WE)
      target_link_libraries(${target} PUBLIC CGAL::Eigen_support)
    endif()
  endif()
endforeach()
