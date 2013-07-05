# This is the CMake script for compiling this folder.

project( Surface_reconstruction_points_3_example )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

# Require packages new or improved since CGAL 3.5 beta 1
include_directories (BEFORE ../../../Installation/include)
include_directories (BEFORE ../../../Point_set_processing_3/include)
include_directories (BEFORE ../../../Surface_mesher/include)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../../../Installation/cmake/modules ${CMAKE_MODULE_PATH})

# Include this package's headers first
include_directories (BEFORE . include ../../include)

# Find CGAL
find_package(CGAL)

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

  # Find Eigen3 (requires 3.1.0 or greater)
  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
  else()
    find_package(TAUCS)
    if(TAUCS_FOUND)
      include( ${TAUCS_USE_FILE} )
    endif(TAUCS_FOUND)
  endif()
  
  if (EIGEN3_FOUND OR TAUCS_FOUND)
    # Executables that require Eigen 3 or TAUCS
    create_single_source_cgal_program( "poisson_reconstruction_example.cpp" )
    create_single_source_cgal_program( "poisson_reconstruction.cpp" )
  else()
    message(STATUS "NOTICE: The examples needs Eigen 3.1 (or greater) or the TAUCS library and will not be compiled.")
  endif()

  
else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

