# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(Arrangement_on_surface_2_Tests)

enable_testing()

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})
# Since CMake-2.8.12: New CMake script, that defines the targets and
# the CTest test cases.
include(${CMAKE_CURRENT_SOURCE_DIR}/cgal_test.cmake)
