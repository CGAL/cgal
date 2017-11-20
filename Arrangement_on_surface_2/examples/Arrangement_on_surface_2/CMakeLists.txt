# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Arrangement_on_surface_2_Examples )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

  if (CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
  endif()
else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
