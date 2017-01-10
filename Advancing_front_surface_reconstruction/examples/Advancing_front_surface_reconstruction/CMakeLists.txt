# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Advancing_front_surface_reconstruction_Examples )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  if (MSVC AND ( CMAKE_SIZEOF_VOID_P EQUAL 4 ) )
     SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244")
  endif()

  include_directories (BEFORE "../../include")

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

