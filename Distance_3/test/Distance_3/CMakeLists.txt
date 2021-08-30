# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Distance_3_Tests)

find_package(CGAL REQUIRED)

if(MSVC)
  # Turn off a VC++ warning on a potential division by zero
  # in Cartesian_kernel/include/CGAL/constructions/kernel_ftC3.h
  # where CGAL_assume() does not help
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4723")
endif()

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()
