# This is the CMake script for compiling this folder.

project( Poisson_surface_reconstruction_3_Examples )

cmake_minimum_required(VERSION 3.1)

# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Allow Windows 32bit applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

    # Print new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # Find Eigen3 (requires 3.1.0 or greater)
  find_package(Eigen3 3.1.0)
  if (EIGEN3_FOUND)
    include( ${EIGEN3_USE_FILE} )
  endif()
  
  if (EIGEN3_FOUND)
    # Executables that require Eigen 3
    create_single_source_cgal_program( "poisson_reconstruction_example.cpp" )
    create_single_source_cgal_program( "poisson_reconstruction.cpp" )
    create_single_source_cgal_program( "poisson_reconstruction_function.cpp" )
  else()
    message(STATUS "NOTICE: The examples need Eigen 3.1 (or greater) will not be compiled.")
  endif()

  
else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

