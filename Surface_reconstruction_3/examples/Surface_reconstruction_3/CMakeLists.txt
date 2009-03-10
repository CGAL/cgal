# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Surface_reconstruction_3_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
# Find CGAL
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # Optimization for applications dealing with large data
  if (MSVC)
    # Use /EHa option to catch stack overflows.
    # Note: TAUCS needs a stack >= 2MB. CGAL default is 10MB.
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHa")

    # Allow Windows applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Turn off stupid VC++ warnings
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4311 /wd4800 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")
    
    # Speed up debugger
    SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_HAS_ITERATOR_DEBUGGING=0")

    # Print new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  include_directories (BEFORE ../../include)
  include_directories (BEFORE include)

  # Temporary debugging stuff
  #ADD_DEFINITIONS( "-DDEBUG_TRACE" )

  # Executables that do *not* require BLAS, LAPACK nor TAUCS
  create_single_source_cgal_program( "APSS_reconstruction.cpp" )
  create_single_source_cgal_program( "average_spacing_example.cpp" )
  create_single_source_cgal_program( "outlier_removal_example.cpp" )
  create_single_source_cgal_program( "pca_normal_estimation_example.cpp" )
  create_single_source_cgal_program( "random_simplification_example.cpp" )
  create_single_source_cgal_program( "surface_reconstruction_read_write_xyz_example.cpp" )

  # Executables that require BLAS and LAPACK
  find_package(LAPACK)
  if(LAPACK_FOUND)

    include( ${LAPACK_USE_FILE} )

    create_single_source_cgal_program( "jet_smoothing_example.cpp" )
    create_single_source_cgal_program( "normal_estimation.cpp" )

  else(LAPACK_FOUND)

    message(STATUS "NOTICE: Some of the executables in this directory need the LAPACK library and will not be compiled.")

  endif(LAPACK_FOUND)

  # Executables that require BLAS, LAPACK and TAUCS 
  find_package(TAUCS)
  if(TAUCS_FOUND)

    include( ${TAUCS_USE_FILE} )

    create_single_source_cgal_program( "poisson_reconstruction.cpp" )
    create_single_source_cgal_program( "poisson_reconstruction_example.cpp" )

  else(TAUCS_FOUND)

    message(STATUS "NOTICE: Some of the executables in this directory need the TAUCS library and will not be compiled.")

  endif(TAUCS_FOUND)

else()
  
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
  
endif()

