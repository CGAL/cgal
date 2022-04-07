# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.20)
project(Polygon_mesh_processing)

# CGAL and its components
find_package(CGAL REQUIRED)

# Boost and its components
find_package(Boost REQUIRED)

if(NOT Boost_FOUND)

  message(
    STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package
find_package(Eigen3 REQUIRED 3.2.0) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

set(FAST_ENVELOPE_BUILD_DIR "" CACHE PATH "Path to fast-evelope build directory")
if (FAST_ENVELOPE_BUILD_DIR)
  message(STATUS "Using ${FAST_ENVELOPE_BUILD_DIR} as build directory of fast-evelope")

  include_directories("${FAST_ENVELOPE_BUILD_DIR}/include")
  link_directories ( "${FAST_ENVELOPE_BUILD_DIR}/lib"  "${FAST_ENVELOPE_BUILD_DIR}"  "${FAST_ENVELOPE_BUILD_DIR}/tbb"   )

  find_package(OpenMP)
  if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()

  create_single_source_cgal_program( "fastE.cpp" )
  target_link_libraries( fastE PUBLIC CGAL::Eigen3_support)
  target_link_libraries( fastE PUBLIC FastEnvelope IndirectPredicates geogram)

else()
  message(STATUS "Cmake variable FAST_ENVELOPE_BUILD_DIR is not defined fastE will not be built")
endif()

create_single_source_cgal_program( "fast.cpp" )
create_single_source_cgal_program("polygon_mesh_slicer.cpp")


if(TARGET CGAL::Eigen3_support)
  target_link_libraries(polygon_mesh_slicer PUBLIC CGAL::Eigen3_support)
endif()
