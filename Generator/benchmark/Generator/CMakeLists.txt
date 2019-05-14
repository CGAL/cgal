# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Generator_example )

cmake_minimum_required(VERSION 2.8.11)


find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  find_package(Boost REQUIRED program_options)
  include_directories(${Boost_PROGRAM_OPTIONS_INCLUDE_DIR})
  add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
  list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_PROGRAM_OPTIONS_LIBRARY})
 
  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../../include)

  create_single_source_cgal_program( "random_grid.cpp" )
  create_single_source_cgal_program( "random_disc_2.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

