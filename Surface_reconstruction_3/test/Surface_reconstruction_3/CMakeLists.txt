# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_reconstruction_3_test )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE include)
  include_directories (BEFORE ../../include)

    create_single_source_cgal_program( "analysis_test.cpp" )
    create_single_source_cgal_program( "APSS_reconstruction_test.cpp" )
    create_single_source_cgal_program( "outlier_removal_test.cpp" )
  # Executables that require LAPACK
  find_package(LAPACK QUIET)
  if(LAPACK_FOUND)

    include( ${LAPACK_USE_FILE} )

    create_single_source_cgal_program( "normal_estimation_test.cpp" )
    create_single_source_cgal_program( "smoothing_test.cpp" )

  else(LAPACK_FOUND)

    message(STATUS "Notice: Some of the executables in this directory need the LAPACK library and will not be compiled.")

  endif(LAPACK_FOUND)

  # Executables that require TAUCS
  find_package(TAUCS QUIET)
  if(TAUCS_FOUND)

    include( ${TAUCS_USE_FILE} )

    create_single_source_cgal_program( "poisson_reconstruction_test.cpp" )

  else(TAUCS_FOUND)

    message(STATUS "Notice: Some of the executables in this directory need the TAUCS library and will not be compiled.")

  endif(TAUCS_FOUND)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

