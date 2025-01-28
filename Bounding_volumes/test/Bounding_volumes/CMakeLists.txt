# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Bounding_volumes_Tests)

find_package(CGAL REQUIRED COMPONENTS Core)

# Use Eigen
find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  if(NOT (${cppfile} STREQUAL "Approximate_min_ellipsoid_d.cpp")
     OR TARGET CGAL::Eigen3_support)
    get_filename_component(target ${cppfile} NAME_WE)
    create_single_source_cgal_program(${cppfile})
    if(TARGET CGAL::Eigen3_support)
      target_link_libraries(${target} PRIVATE CGAL::CGAL CGAL::Eigen3_support)
    else()
      target_link_libraries(${target} PRIVATE CGAL::CGAL)
    endif()
  endif()
endforeach()
