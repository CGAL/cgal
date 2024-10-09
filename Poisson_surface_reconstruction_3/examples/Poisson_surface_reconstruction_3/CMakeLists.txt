# This is the CMake script for compiling this folder.

cmake_minimum_required(VERSION 3.12...3.29)
project(Poisson_surface_reconstruction_3_Examples)

# Find CGAL
find_package(CGAL REQUIRED)

# VisualC++ optimization for applications dealing with large data
if(MSVC)
  # Allow Windows 32bit applications to use up to 3GB of RAM
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")

  # Print new compilation options
  message(STATUS "USING DEBUG CXXFLAGS   = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}'")
  message(STATUS "USING DEBUG EXEFLAGS   = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}'")
  message(STATUS "USING RELEASE CXXFLAGS = '${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}'")
  message(STATUS "USING RELEASE EXEFLAGS = '${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}'")
endif()

# Activate Parallelism in Mesh_3
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)

# Find Eigen3 (requires 3.1.0 or greater)
find_package(Eigen3 3.1.0 QUIET)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  # Executables that require Eigen 3
  create_single_source_cgal_program("poisson_reconstruction_example.cpp")
  target_link_libraries(poisson_reconstruction_example PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("poisson_and_parallel_mesh_3.cpp")
  target_link_libraries(poisson_and_parallel_mesh_3 PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("poisson_reconstruction.cpp")
  target_link_libraries(poisson_reconstruction PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("poisson_reconstruction_function.cpp")
  target_link_libraries(poisson_reconstruction_function PUBLIC CGAL::Eigen3_support)
  create_single_source_cgal_program("tutorial_example.cpp")
  target_link_libraries(tutorial_example PUBLIC CGAL::Eigen3_support)
else()
  message("NOTICE: The examples require Eigen 3.1 (or greater) will not be compiled.")
endif()

# Use TBB for parallelism
if(CGAL_ACTIVATE_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
  if(TARGET CGAL::TBB_support)
    add_definitions(-DCGAL_CONCURRENT_MESH_3)
    target_link_libraries(poisson_and_parallel_mesh_3 PUBLIC CGAL::TBB_support)
  endif()#Find TBB
endif()#Parallelism in Mesh_3
