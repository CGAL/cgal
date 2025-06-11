# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Arrangement_on_surface_2_Tests)

enable_testing()

find_package(CGAL REQUIRED COMPONENTS Core Qt6)

include(${CMAKE_CURRENT_SOURCE_DIR}/cgal_test.cmake)
