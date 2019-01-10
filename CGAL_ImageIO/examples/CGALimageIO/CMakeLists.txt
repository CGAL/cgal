# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( CGALimageIO_Examples ) 

cmake_minimum_required(VERSION 3.1)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(CGAL QUIET COMPONENTS ImageIO )

if(CGAL_ImageIO_FOUND)
  include(${CGAL_USE_FILE})

  create_single_source_cgal_program( "convert_raw_image_to_inr.cpp" )
  create_single_source_cgal_program( "test_imageio.cpp" )
else()
  message(STATUS "NOTICE: This demo needs the CGAL ImageIO library, and will not be compiled.")
endif()
