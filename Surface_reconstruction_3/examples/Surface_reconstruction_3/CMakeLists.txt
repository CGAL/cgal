# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_reconstruction_3_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CreateSingleSourceCGALProgram )

  include_directories (BEFORE include)
  include_directories (BEFORE ../../include)

  create_single_source_cgal_program( "APSS_reconstruction.cpp" )
  create_single_source_cgal_program( "average_spacing_example.cpp" )
  create_single_source_cgal_program( "outlier_removal_example.cpp" )
  create_single_source_cgal_program( "pca_normal_estimation_example.cpp" )
  create_single_source_cgal_program( "random_simplification_example.cpp" )
  create_single_source_cgal_program( "surface_reconstruction_read_write_xyz_example.cpp" )

  # Executables that require LAPACK
  find_package(LAPACK QUIET)
  if(LAPACK_FOUND)

    include( ${LAPACK_USE_FILE} )

    create_single_source_cgal_program( "jet_smoothing_example.cpp" )
    create_single_source_cgal_program( "normal_estimation.cpp" )

  else(LAPACK_FOUND)

    message(STATUS "Notice: Some of the executables in this directory need the LAPACK library and will not be compiled.")

  endif(LAPACK_FOUND)

  # Executables that require TAUCS
  find_package(TAUCS QUIET)
  if(TAUCS_FOUND)

    include( ${TAUCS_USE_FILE} )

    create_single_source_cgal_program( "poisson_reconstruction.cpp" )
    create_single_source_cgal_program( "poisson_reconstruction_example.cpp" )

  else(TAUCS_FOUND)

    message(STATUS "Notice: Some of the executables in this directory need the TAUCS library and will not be compiled.")

  endif(TAUCS_FOUND)

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

