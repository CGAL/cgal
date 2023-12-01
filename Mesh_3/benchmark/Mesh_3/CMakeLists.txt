# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Mesh_3_benchmark)

include_directories(../../../Triangulation_3/include)
include_directories(../../../STL_Extension/include)
include_directories(../../../AABB_tree/include)
add_definitions(-DCGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
                -DCGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS)

if(MESH_3_VERBOSE)
  add_definitions(-DCGAL_MESH_3_VERBOSE)
endif()

find_package(CGAL REQUIRED COMPONENTS ImageIO)

# Activate concurrency ? (turned OFF by default)
option(CGAL_ACTIVATE_CONCURRENT_MESH_3 "Activate parallelism in Mesh_3" OFF)

# And add -DCGAL_CONCURRENT_MESH_3 if that option is ON
if(CGAL_ACTIVATE_CONCURRENT_MESH_3)
  add_definitions(-DCGAL_CONCURRENT_MESH_3)
  find_package(TBB REQUIRED)
  include(CGAL_TBB_support)
else()
  option(LINK_WITH_TBB
         "Link with TBB anyway so we can use TBB timers for profiling" ON)
  if(LINK_WITH_TBB)
    find_package(TBB)
    include(CGAL_TBB_support)
  endif(LINK_WITH_TBB)
endif()

# Compilable benchmark
set(BENCHMARK_SOURCE_FILES "concurrency.cpp")
add_msvc_precompiled_header("StdAfx.h" "StdAfx.cpp" BENCHMARK_SOURCE_FILES)
create_single_source_cgal_program(${BENCHMARK_SOURCE_FILES})
if(TARGET CGAL::TBB_support)
  target_link_libraries(concurrency PUBLIC CGAL::TBB_support)
endif()

# Link with Boost.ProgramOptions (optional)
find_package(Boost QUIET COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
  if(TARGET Boost::program_options)
    target_link_libraries(concurrency PRIVATE Boost::program_options)
  else()
    target_link_libraries(concurrency PRIVATE ${Boost_PROGRAM_OPTIONS_LIBRARY})
  endif()
endif()
