# This is the CMake script for compiling this folder.

project( Point_set_processing_3_test )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

# Require packages new or improved since CGAL 3.4
include_directories (BEFORE ../../../Installation/include/)

# Include this package's headers first
include_directories (BEFORE . include ../../include)

# Find CGAL
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Use /FR to turn on IntelliSense
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /FR")

    # Allow Windows applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Turn off stupid VC++ warnings
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4267 /wd4311 /wd4800 /wd4503 /wd4244 /wd4345 /wd4996 /wd4396 /wd4018")

    # Prints new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # Temporary debugging stuff
  ADD_DEFINITIONS( "-DDEBUG_TRACE" )

  find_package( TBB QUIET )

  if( TBB_FOUND )
    include(${TBB_USE_FILE})
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
  endif()
  
  # Executables that do *not* require LAPACK
  create_single_source_cgal_program( "read_test.cpp" )
  create_single_source_cgal_program( "read_test_with_different_pmaps.cpp" )
  create_single_source_cgal_program( "analysis_test.cpp" )
  create_single_source_cgal_program( "remove_outliers_test.cpp" )
  create_single_source_cgal_program( "wlop_simplify_and_regularize_test.cpp" )
  create_single_source_cgal_program( "bilateral_smoothing_test.cpp" )
  create_single_source_cgal_program( "edge_aware_upsample_test.cpp" )

  # Use Eigen or BLAS and LAPACK (optional)
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if (NOT EIGEN3_FOUND)
    find_package(LAPACK)
    if(LAPACK_FOUND)
      include( ${LAPACK_USE_FILE} )
    endif(LAPACK_FOUND)
  else()
    include( ${EIGEN3_USE_FILE} )
  endif()
  
  
  if(EIGEN3_FOUND OR LAPACK_FOUND)
    # Executables that require Eigen or BLAS and LAPACK
    create_single_source_cgal_program( "normal_estimation_test.cpp" )
    create_single_source_cgal_program( "smoothing_test.cpp" )

  else(EIGEN3_FOUND OR LAPACK_FOUND)

    message(STATUS "NOTICE: This program requires either Eigen 3.1 (or greater) or LAPACK, and will not be compiled.")

  endif(EIGEN3_FOUND OR LAPACK_FOUND)

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

