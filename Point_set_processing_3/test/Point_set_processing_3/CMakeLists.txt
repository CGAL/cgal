# This is the CMake script for compiling this folder.

project( Point_set_processing_3_test )

cmake_minimum_required(VERSION 2.8.11)


# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    # Allow Windows 32bit applications to use up to 3GB of RAM
    SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

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
    create_single_source_cgal_program( "hierarchy_simplification_test.cpp" )
    create_single_source_cgal_program( "smoothing_test.cpp" )
    create_single_source_cgal_program( "vcm_plane_test.cpp" )
    create_single_source_cgal_program( "vcm_all_test.cpp" )
    create_single_source_cgal_program( "jet_pointer_as_property_map.cpp" )

  else(EIGEN3_FOUND OR LAPACK_FOUND)

    message(STATUS "NOTICE: This program requires either Eigen 3.1 (or greater) or LAPACK, and will not be compiled.")

  endif(EIGEN3_FOUND OR LAPACK_FOUND)

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

