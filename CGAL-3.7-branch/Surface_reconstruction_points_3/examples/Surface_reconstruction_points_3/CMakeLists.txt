# This is the CMake script for compiling this folder.

project( Surface_reconstruction_points_3_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )
endif(COMMAND cmake_policy)

# Require packages new or improved since CGAL 3.5 beta 1
include_directories (BEFORE ../../../Installation/include)
include_directories (BEFORE ../../../Point_set_processing_3/include)
include_directories (BEFORE ../../../Surface_mesher/include)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../../Installation/cmake/modules ${CMAKE_MODULE_PATH})

# TEMPORARY: require experimental-packages/ARAP_Parameterization for MKL support
include_directories (BEFORE ../../../../experimental-packages/ARAP_Parameterization/include)

# Include this package's headers first
include_directories (BEFORE . include ../../include)

# Find CGAL
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Use /EHa option to catch stack overflows.
    # Note: TAUCS needs a stack >= 2MB. CGAL default is 10MB.
    string(REGEX REPLACE "/EH[asc]*" "/EHa" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

    # Use /FR to turn on IntelliSense
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FR")

    # Allow Windows applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Turn off stupid VC++ warnings
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4311 /wd4800 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")

    # Print new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # Temporary debugging stuff
  ADD_DEFINITIONS( "-DDEBUG_TRACE" ) # turn on traces

  # Link with BLAS, LAPACK and TAUCS (optional)
  find_package(TAUCS)
  if(TAUCS_FOUND)

    include( ${TAUCS_USE_FILE} )

    # Executables that require TAUCS but not MKL
    create_single_source_cgal_program( "poisson_reconstruction_example.cpp" )

    # Link with TAUCS and MKL (optional)
    find_package(MKL)
    if(MKL_FOUND)

      include( ${MKL_USE_FILE} )

      # Executables that require TAUCS and MKL
      create_single_source_cgal_program( "poisson_reconstruction.cpp" )

    else(MKL_FOUND)

      message(STATUS "NOTICE: Some of the executables in this directory need the MKL library and will not be compiled.")

    endif(MKL_FOUND)

  else(TAUCS_FOUND)

    message(STATUS "NOTICE: Some of the executables in this directory need the TAUCS library and will not be compiled.")

  endif(TAUCS_FOUND)

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

