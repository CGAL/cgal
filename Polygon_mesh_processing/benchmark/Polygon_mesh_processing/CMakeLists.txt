# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)
project(Polygon_mesh_processing)

# CGAL and its components
find_package(CGAL REQUIRED)

find_package(Eigen3 3.2.0 QUIET) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: Benchmarks require Eigen 3.2 (or greater), and will not be compiled")
  return()
endif()

set(FAST_ENVELOPE_BUILD_DIR "" CACHE PATH "Path to fast-evelope build directory")
if (FAST_ENVELOPE_BUILD_DIR)
  message(STATUS "Using ${FAST_ENVELOPE_BUILD_DIR} as build directory of fast-evelope")

  include_directories("${FAST_ENVELOPE_BUILD_DIR}/include")
  link_directories ( "${FAST_ENVELOPE_BUILD_DIR}/lib" "${FAST_ENVELOPE_BUILD_DIR}" "${FAST_ENVELOPE_BUILD_DIR}/tbb")

  find_package(OpenMP)
  if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()

  create_single_source_cgal_program("fastE.cpp")
  target_link_libraries( fastE PUBLIC CGAL::Eigen3_support)
  target_link_libraries( fastE PUBLIC FastEnvelope IndirectPredicates geogram)

else()
  message(STATUS "CMake variable FAST_ENVELOPE_BUILD_DIR is not defined; benchmark 'fastE' will not be built")
endif()

create_single_source_cgal_program("fast.cpp")

create_single_source_cgal_program("polygon_mesh_slicer.cpp")
target_link_libraries(polygon_mesh_slicer PUBLIC CGAL::Eigen3_support)

