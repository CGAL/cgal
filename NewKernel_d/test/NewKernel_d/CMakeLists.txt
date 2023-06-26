# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.
# Then modified by hand to add Eigen3.

cmake_minimum_required(VERSION 3.1...3.23)
project(NewKernel_d_Tests)

if(CMAKE_COMPILER_IS_GNUCCX AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
  message("NOTICE: this directory requires a version of gcc >= 4.4, and will not be compiled.")
  return()
endif()

find_package(CGAL REQUIRED)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  file(
    GLOB cppfiles
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    get_filename_component(target ${cppfile} NAME_WE)
    create_single_source_cgal_program("${cppfile}")
    target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
  endforeach()
else()
  message("NOTICE: These programs require the Eigen3 library, and will not be compiled.")
endif()
