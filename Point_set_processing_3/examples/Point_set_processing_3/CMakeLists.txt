# This is the CMake script for compiling this folder.

project( Point_set_processing_3_Examples )

cmake_minimum_required(VERSION 2.8.11)


# Find CGAL
find_package(CGAL QUIET)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )

  find_package(Boost QUIET)

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

  # Activate concurrency?
  option(ACTIVATE_CONCURRENT_PSP3
         "Enable concurrency"
         OFF)
     
  if( ACTIVATE_CONCURRENT_PSP3 OR ENV{ACTIVATE_CONCURRENT_PSP3} )
    find_package( TBB REQUIRED )
    if( TBB_FOUND )
      include(${TBB_USE_FILE})
      list(APPEND CGAL_3RD_PARTY_LIBRARIES ${TBB_LIBRARIES})
    endif()
  endif()

  # Executables that do *not* require EIGEN or LAPACK
  create_single_source_cgal_program( "average_spacing_example.cpp" )
  create_single_source_cgal_program( "bilateral_smooth_point_set_example.cpp" )
  create_single_source_cgal_program( "grid_simplification_example.cpp" )
  create_single_source_cgal_program( "grid_simplify_indices.cpp" )
  create_single_source_cgal_program( "hierarchy_simplification_example.cpp" )
  create_single_source_cgal_program( "normals_example.cpp" )
  create_single_source_cgal_program( "property_map.cpp" )
  create_single_source_cgal_program( "random_simplification_example.cpp" )
  create_single_source_cgal_program( "read_write_xyz_point_set_example.cpp" )
  create_single_source_cgal_program( "remove_outliers_example.cpp" )
  create_single_source_cgal_program( "scale_estimation_example.cpp" )
  create_single_source_cgal_program( "scale_estimation_2d_example.cpp" )
  create_single_source_cgal_program( "wlop_simplify_and_regularize_point_set_example.cpp" )
  create_single_source_cgal_program( "edge_aware_upsample_point_set_example.cpp" )
  create_single_source_cgal_program( "structuring_example.cpp" )

  set(needed_cxx_features cxx_rvalue_references cxx_variadic_templates)    
  create_single_source_cgal_program( "read_ply_points_with_colors_example.cpp" CXX_FEATURES ${needed_cxx_features} )
  create_single_source_cgal_program( "write_ply_points_example.cpp" CXX_FEATURES ${needed_cxx_features} )
    
  find_package(LASLIB)
  if (LASLIB_FOUND)
    include(${LASLIB_USE_FILE})
    include_directories(${LASLIB_INCLUDE_DIR})
    include_directories(${LASZIP_INCLUDE_DIR})
    create_single_source_cgal_program( "read_las_example.cpp" CXX_FEATURES ${needed_cxx_features} )
    target_link_libraries(read_las_example ${LASLIB_LIBRARIES})
  else()
    message(STATUS "NOTICE : the LAS reader test requires LASlib and will not be compiled.")
  endif()
  
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
    create_single_source_cgal_program( "jet_smoothing_example.cpp" )
    create_single_source_cgal_program( "normal_estimation.cpp" )
    create_single_source_cgal_program( "edges_example.cpp" )
  else(EIGEN3_FOUND OR LAPACK_FOUND)

    message(STATUS "NOTICE: Some of the executables in this directory need either Eigen 3.1 (or greater) or LAPACK, and will not be compiled.")

  endif(EIGEN3_FOUND OR LAPACK_FOUND)


else()
    message(STATUS "NOTICE: This program requires the CGAL library, and will not be compiled.")
endif()
