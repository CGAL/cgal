# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Jet_fitting_3_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  find_package(LAPACK QUIET)
  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )

    # Link with Boost.ProgramOptions (optional)
    find_package(Boost QUIET COMPONENTS program_options)
    if(Boost_PROGRAM_OPTIONS_FOUND)
      add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
      link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
    endif()

    create_single_source_cgal_program( "Mesh_estimation.cpp" )
    create_single_source_cgal_program( "Single_estimation.cpp" )

  else(LAPACK_FOUND)

    message(STATUS "NOTICE: This program requires the LAPACK API, and will not be compiled.")

  endif(LAPACK_FOUND)

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

