# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( CGALImageIO_example ) 

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS ImageIO )

include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )

if(CGAL_ImageIO_FOUND)
  create_single_source_cgal_program( "convert_raw_image_to_inr.cpp" )
  create_single_source_cgal_program( "test_imageio.cpp" )
  if(CGAL_ImageIO_USE_ZLIB)
    target_link_libraries(convert_raw_image_to_inr ${ZLIB_LIBRARY})
    target_link_libraries(test_imageio ${ZLIB_LIBRARY})
  endif()
else()
  message(STATUS "NOTICE: This demo needs the CGAL ImageIO library, and will not be compiled.")
endif()
