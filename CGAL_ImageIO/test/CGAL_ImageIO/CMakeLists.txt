# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.20)
project(CGAL_ImageIO_Tests)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(CGAL REQUIRED COMPONENTS ImageIO)

include(${CGAL_USE_FILE})

if(WITH_CGAL_ImageIO)
  create_single_source_cgal_program("test_trilinear_interpolation.cpp")
else()
  message(
    STATUS
      "NOTICE: Some tests require the CGAL_ImageIO library, and will not be compiled."
  )
endif()
