# This is the CMake script for compiling a CGAL application.
project( Ridges_3_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  find_package(LAPACK QUIET)

  if(LAPACK_FOUND)
    include( ${LAPACK_USE_FILE} )

    find_package( Boost QUIET COMPONENTS program_options )
    if( Boost_PROGRAM_OPTIONS_FOUND )
      add_definitions( -DCGAL_USE_BOOST_PROGRAM_OPTIONS )
      link_libraries( ${Boost_PROGRAM_OPTIONS_LIBRARY} )
    endif( Boost_PROGRAM_OPTIONS_FOUND )

    add_executable(Compute_Ridges_Umbilics Compute_Ridges_Umbilics.cpp PolyhedralSurf.cpp)
    target_link_libraries(Compute_Ridges_Umbilics ${CGAL_LIBRARIES})

  else(LAPACK_FOUND)

    message(STATUS "NOTICE: This program needs the LAPACK API and will not be compiled.")

  endif(LAPACK_FOUND)

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()
