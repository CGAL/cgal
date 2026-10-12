# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.12...3.31)
project(Point_set_processing_3_Examples)

# Find CGAL
find_package(CGAL REQUIRED)

# VisualC++ optimization for applications dealing with large data
if(MSVC)
  if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # Allow Windows 32bit applications to use up to 3GB of RAM
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
  endif()
  # Prints new compilation options
  message(STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'")
  message(STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'")
  message(STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'")
  message(STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'")
endif()

# Activate concurrency?
option(CGAL_ACTIVATE_CONCURRENT_PSP3 "Enable concurrency" OFF)

if(CGAL_ACTIVATE_CONCURRENT_PSP3 OR "$ENV{CGAL_ACTIVATE_CONCURRENT_PSP3}")
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
endif()

# Executables that do *not* require EIGEN
foreach(
  target
  average_spacing_example
  bilateral_smooth_point_set_example
  grid_simplification_example
  grid_simplify_indices
  property_map
  random_simplification_example
  read_write_xyz_point_set_example
  remove_outliers_example
  wlop_simplify_and_regularize_point_set_example
  edge_aware_upsample_point_set_example
  structuring_example
  read_ply_points_with_colors_example
  write_ply_points_example
  poisson_eliminate_example
  poisson_eliminate_from_mesh_example)
  create_single_source_cgal_program("${target}.cpp")
  target_link_libraries(${target} PRIVATE CGAL::CGAL)
endforeach()

#disable if MSVC 2017
if(NOT MSVC_VERSION OR MSVC_VERSION GREATER_EQUAL 1919 OR MSVC_VERSION LESS 1910)
  find_package(LASLIB QUIET)
  include(CGAL_LASLIB_support)
  if (TARGET CGAL::LASLIB_support)
    create_single_source_cgal_program("read_las_example.cpp")
    create_single_source_cgal_program("write_las_example.cpp")
    target_link_libraries(read_las_example PRIVATE CGAL::LASLIB_support)
    target_link_libraries(write_las_example PRIVATE CGAL::LASLIB_support)

    target_compile_definitions(read_las_example PRIVATE -D_CRT_SECURE_NO_DEPRECATE -D_CRT_SECURE_NO_WARNINGS)
    target_compile_definitions(write_las_example PRIVATE -D_CRT_SECURE_NO_DEPRECATE -D_CRT_SECURE_NO_WARNINGS)

  else()
    message(STATUS "NOTICE : the LAS reader test requires LASlib and will not be compiled.")
  endif()
else()
  message(STATUS "NOTICE: the LAS reader does not work with your version of Visual Studio 2017.")
endif()

# Use Eigen
find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)

  # Executables that require Eigen
  foreach(
    target
    callback_example
    clustering_example
    edges_example
    hierarchy_simplification_example
    jet_smoothing_example
    normal_estimation
    normals_example
    scale_estimation_example
    scale_estimation_2d_example)
    create_single_source_cgal_program("${target}.cpp")
    target_link_libraries(${target} PRIVATE CGAL::Eigen3_support)
  endforeach()

  if (TARGET CGAL::TBB_support)
    foreach(
      target
      callback_example
      clustering_example
      jet_smoothing_example
      normal_estimation
      normals_example
      scale_estimation_example)
      target_link_libraries(${target} PRIVATE CGAL::TBB_support)
    endforeach()
  endif()

  if (TARGET CGAL::Eigen3_support AND TARGET CGAL::LASLIB_support)
    create_single_source_cgal_program( "orient_scanlines_example.cpp" )
    target_link_libraries(orient_scanlines_example PRIVATE
      CGAL::Eigen3_support CGAL::LASLIB_support)
    if (TARGET CGAL::TBB_support)
      target_link_libraries(orient_scanlines_example PRIVATE CGAL::TBB_support)
    endif()
  endif()

  # Executables that require libpointmatcher
  find_package(libpointmatcher QUIET)
  include(CGAL_pointmatcher_support)
  if(TARGET CGAL::pointmatcher_support)
    create_single_source_cgal_program("registration_with_pointmatcher.cpp")
    target_link_libraries(registration_with_pointmatcher PRIVATE CGAL::Eigen3_support
                          CGAL::pointmatcher_support)
  else()
    message(STATUS "NOTICE: registration with pointmatcher requires libpointmatcher and will not be compiled.")
  endif()

  # Executables that require OpenGR
  find_package(OpenGR QUIET)
  include(CGAL_OpenGR_support)

  if(TARGET CGAL::OpenGR_support)
    create_single_source_cgal_program("registration_with_OpenGR.cpp")
    target_link_libraries(registration_with_OpenGR PRIVATE CGAL::Eigen3_support
                          CGAL::OpenGR_support)
  else()
    message(STATUS "NOTICE: registration_with_OpenGR requires OpenGR, and will not be compiled.")
  endif()

  # Executables that require both libpointmatcher and OpenGR
  if (TARGET CGAL::pointmatcher_support AND
      TARGET CGAL::OpenGR_support)
    create_single_source_cgal_program("registration_with_opengr_pointmatcher_pipeline.cpp")
    target_link_libraries(
      registration_with_opengr_pointmatcher_pipeline PRIVATE CGAL::Eigen3_support
      CGAL::pointmatcher_support CGAL::OpenGR_support)
  else()
    message(STATUS "NOTICE: registration with OpenGR and pointmatcher requires both libpointmatcher and OpenGR, and will not be compiled.")
  endif()

else()
  message(STATUS "NOTICE: Some of the executables in this directory require Eigen 3.1 (or greater), and will not be compiled.")
endif()

if(CGAL_ACTIVATE_CONCURRENT_MESH_3 AND TARGET CGAL::TBB_support)
  foreach(
    target
    average_spacing_example
    bilateral_smooth_point_set_example
    edge_aware_upsample_point_set_example
    remove_outliers_example
    wlop_simplify_and_regularize_point_set_example)
    target_link_libraries(${target} PRIVATE CGAL::TBB_support)
  endforeach()
endif()
