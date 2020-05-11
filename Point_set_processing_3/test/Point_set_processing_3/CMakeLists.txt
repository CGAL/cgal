# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.1...3.15)
project( Point_set_processing_3_Tests )



# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  # VisualC++ optimization for applications dealing with large data
  if (MSVC)
    if ( CMAKE_SIZEOF_VOID_P EQUAL 4 )
      # Allow Windows 32bit applications to use up to 3GB of RAM
      SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
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

  # Executables that do *not* require Eigen
  create_single_source_cgal_program( "read_test.cpp" )
  create_single_source_cgal_program( "test_read_write_point_set.cpp" )
  create_single_source_cgal_program( "read_test_with_different_pmaps.cpp" )
  create_single_source_cgal_program( "analysis_test.cpp" )
  create_single_source_cgal_program( "remove_outliers_test.cpp" )
  create_single_source_cgal_program( "wlop_simplify_and_regularize_test.cpp" )
  create_single_source_cgal_program( "bilateral_smoothing_test.cpp" )
  create_single_source_cgal_program( "edge_aware_upsample_test.cpp" )
  create_single_source_cgal_program( "structuring_test.cpp" )

  #Use LAS
  find_package(LASLIB)
  if (LASLIB_FOUND)
    CGAL_target_use_LASLIB(test_read_write_point_set)
  else()
    message(STATUS "NOTICE : the LAS reader test requires LASlib and will not be compiled.")
  endif()

  # Use Eigen
  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if (EIGEN3_FOUND)
    # Executables that require Eigen
    create_single_source_cgal_program( "normal_estimation_test.cpp" )
    CGAL_target_use_Eigen(normal_estimation_test)

    create_single_source_cgal_program( "hierarchy_simplification_test.cpp" )
    CGAL_target_use_Eigen(hierarchy_simplification_test)

    create_single_source_cgal_program( "smoothing_test.cpp" )
    CGAL_target_use_Eigen(smoothing_test)

    create_single_source_cgal_program( "vcm_plane_test.cpp" )
    CGAL_target_use_Eigen(vcm_plane_test)

    create_single_source_cgal_program( "vcm_all_test.cpp" )
    CGAL_target_use_Eigen(vcm_all_test)

    create_single_source_cgal_program( "jet_pointer_as_property_map.cpp" )
    CGAL_target_use_Eigen(jet_pointer_as_property_map)
  else()
    message(STATUS "NOTICE: This program requires Eigen 3.1 (or greater) and will not be compiled.")
  endif()

  if (TBB_FOUND)
    foreach(target
        analysis_test
        smoothing_test
        bilateral_smoothing_test
        wlop_simplify_and_regularize_test
        edge_aware_upsample_test
        normal_estimation_test)
      if(TARGET ${target})
        CGAL_target_use_TBB(${target})
      endif()
    endforeach()
  endif()

else()

    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")

endif()

