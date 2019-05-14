# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Cone_spanners_2_Examples )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core)

if ( CGAL_FOUND AND CGAL_Core_FOUND)

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")

  if (MSVC)
    # Turn off a boost related warning that appears with VC2015
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4172")
  endif()

  # create a target per cppfile
  file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
  foreach(cppfile ${cppfiles})
    create_single_source_cgal_program( "${cppfile}" )
  endforeach()

else()

  message(STATUS "This program requires the CGAL and CGAL_Core libraries, and will not be compiled.")

endif()

