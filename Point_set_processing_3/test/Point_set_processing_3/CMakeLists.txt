# This is the CMake script for compiling this folder.

project( Point_set_processing_3_Tests )

cmake_minimum_required(VERSION 2.8.11)


# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    if ( CMAKE_SIZEOF_VOID_P EQUAL 4 )
      # Allow Windows 32bit applications to use up to 3GB of RAM
      SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4244")
    endif()
    # Prints new compilation options
    message( STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'" )
    message( STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'" )
    message( STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'" )
    message( STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'" )
  endif()

  # Temporary debugging stuff
  ADD_DEFINITIONS( "-DDEBUG_TRACE" )

  find_package( TBB QUIET )

  # Executables that do *not* require LAPACK
  create_single_source_cgal_program( "read_test.cpp" )
  create_single_source_cgal_program( "read_test_with_different_pmaps.cpp" )
  create_single_source_cgal_program( "analysis_test.cpp" )
  create_single_source_cgal_program( "remove_outliers_test.cpp" )
  create_single_source_cgal_program( "wlop_simplify_and_regularize_test.cpp" )
  create_single_source_cgal_program( "bilateral_smoothing_test.cpp" )
  create_single_source_cgal_program( "edge_aware_upsample_test.cpp" )
  create_single_source_cgal_program( "structuring_test.cpp" )

  # Deprecated tests
  create_single_source_cgal_program( "deprecated_read_test.cpp" )
  create_single_source_cgal_program( "deprecated_read_test_with_different_pmaps.cpp" )
  create_single_source_cgal_program( "deprecated_analysis_test.cpp" )
  create_single_source_cgal_program( "deprecated_remove_outliers_test.cpp" )
  create_single_source_cgal_program( "deprecated_wlop_simplify_and_regularize_test.cpp" )
  create_single_source_cgal_program( "deprecated_bilateral_smoothing_test.cpp" )
  create_single_source_cgal_program( "deprecated_edge_aware_upsample_test.cpp" )
  create_single_source_cgal_program( "deprecated_structuring_test.cpp" )

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

    # Deprecated tests
    create_single_source_cgal_program( "deprecated_normal_estimation_test.cpp" )
    create_single_source_cgal_program( "deprecated_hierarchy_simplification_test.cpp" )
    create_single_source_cgal_program( "deprecated_smoothing_test.cpp" )
    create_single_source_cgal_program( "deprecated_vcm_plane_test.cpp" )
    create_single_source_cgal_program( "deprecated_vcm_all_test.cpp" )
    create_single_source_cgal_program( "deprecated_jet_pointer_as_property_map.cpp" )

  else(EIGEN3_FOUND OR LAPACK_FOUND)

    message(STATUS "NOTICE: This program requires either Eigen 3.1 (or greater) or LAPACK, and will not be compiled.")

  endif(EIGEN3_FOUND OR LAPACK_FOUND)

  foreach(target
      analysis_test
      smoothing_test
      bilateral_smoothing_test
      wlop_simplify_and_regularize_test
      edge_aware_upsample_test
      normal_estimation_test)
    if(TBB_FOUND AND TARGET ${target})
      CGAL_target_use_TBB(${target})
    endif()
  endforeach()

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

