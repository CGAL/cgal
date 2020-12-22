# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.1...3.15)
project(Point_set_processing_3_Examples)

# Find CGAL
find_package(CGAL REQUIRED)

find_package(Boost QUIET)

# VisualC++ optimization for applications dealing with large data
if(MSVC)
  # Quit warning in the lasreader
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_SECURE_NO_WARNINGS)

  if(CMAKE_SIZEOF_VOID_P EQUAL 4)
    # Allow Windows 32bit applications to use up to 3GB of RAM
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
  endif()
  # Prints new compilation options
  message(
    STATUS
      "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'")
  message(
    STATUS
      "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'"
  )
  message(
    STATUS
      "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'"
  )
  message(
    STATUS
      "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'"
  )
endif()

# Activate concurrency?
option(CGAL_ACTIVATE_CONCURRENT_PSP3 "Enable concurrency" OFF)

set(CGAL_libs CGAL::CGAL)
if(CGAL_ACTIVATE_CONCURRENT_PSP3 OR ENV{CGAL_ACTIVATE_CONCURRENT_PSP3})
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
  if(TARGET CGAL::TBB_support)
    set(CGAL_TBB_target ${CGAL_libs} CGAL::TBB_support)
  endif()
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
  write_ply_points_example)
  add_executable(${target} "${target}.cpp")
  target_link_libraries(${target} ${CGAL_libs})
endforeach()

find_package(LASLIB)
include(CGAL_LASLIB_support)
if(TARGET CGAL::LASLIB_support)
  add_executable(read_las_example "read_las_example.cpp")
  target_link_libraries(read_las_example ${CGAL_libs} CGAL::LASLIB_support)
else()
  message(
    STATUS
      "NOTICE : the LAS reader test requires LASlib and will not be compiled.")
endif()

# Use Eigen
find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen_support)
if(TARGET CGAL::Eigen_support)
  set(CGAL_libs ${CGAL_libs} CGAL::Eigen_support)

  # Executables that require Eigen
  foreach(
    target
    jet_smoothing_example
    normal_estimation
    clustering_example
    edges_example
    callback_example
    scale_estimation_example
    scale_estimation_2d_example
    hierarchy_simplification_example
    normals_example)
    add_executable(${target} "${target}.cpp")
    target_link_libraries(${target} ${CGAL_libs})
  endforeach()
  if (TARGET CGAL::Eigen_support AND TARGET CGAL::LASLIB_support)
    add_executable( orient_scanlines_example "orient_scanlines_example.cpp" )
    target_link_libraries(orient_scanlines_example ${CGAL_libs}
      CGAL::Eigen_support CGAL::LASLIB_support)
  endif()
  # Executables that require libpointmatcher
  find_package(libpointmatcher QUIET)
  include(CGAL_pointmatcher_support)
  if(TARGET CGAL::pointmatcher_support)
    add_executable(registration_with_pointmatcher
                   "registration_with_pointmatcher.cpp")
    target_link_libraries(registration_with_pointmatcher ${CGAL_libs}
                          CGAL::pointmatcher_support)
  else()
    message(
      STATUS
        "NOTICE : the registration_with_pointmatcher test requires libpointmatcher and will not be compiled."
    )
  endif()
  # Executables that require OpenGR
  find_package(OpenGR QUIET)
  include(CGAL_OpenGR_support)
  if(TARGET CGAL::OpenGR_support)
    add_executable(registration_with_OpenGR "registration_with_OpenGR.cpp")
    target_link_libraries(registration_with_OpenGR ${CGAL_libs}
                          CGAL::OpenGR_support)
  else()
    message(
      STATUS
        "NOTICE : registration_with_OpenGR requires OpenGR, and will not be compiled."
    )
  endif()

  # Executables that require both libpointmatcher and OpenGR
  if (TARGET CGAL::pointmatcher_support AND
      TARGET CGAL::OpenGR_support)
    add_executable(registration_with_opengr_pointmatcher_pipeline
                   "registration_with_opengr_pointmatcher_pipeline.cpp")
    target_link_libraries(
      registration_with_opengr_pointmatcher_pipeline ${CGAL_libs}
      CGAL::pointmatcher_support CGAL::OpenGR_support)
  else()
    message(
      STATUS
        "NOTICE : registration_with_opengr_pointmatcher_pipeline requires libpointmatcher and OpenGR, and will not be compiled."
    )
  endif()

else()
  message(
    STATUS
      "NOTICE: Some of the executables in this directory need Eigen 3.1 (or greater) and will not be compiled."
  )
endif()
