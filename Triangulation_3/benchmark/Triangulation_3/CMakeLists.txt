# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project( Triangulation_3 )


# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package


# Creating entries for all C++ files with "main" routine
# ##########################################################


create_single_source_cgal_program( "incident_edges.cpp" )

create_single_source_cgal_program( "simple_2.cpp" )

create_single_source_cgal_program( "simple.cpp" )

create_single_source_cgal_program( "Triangulation_benchmark_3.cpp" )

create_single_source_cgal_program( "segment_traverser_benchmark.cpp" )

find_package(benchmark)

if(TARGET benchmark::benchmark)
  find_package(TBB REQUIRED)
  include( CGAL_target_use_TBB )

  create_single_source_cgal_program( "DT3_benchmark_with_TBB.cpp" )
  CGAL_target_use_TBB(DT3_benchmark_with_TBB)
  target_link_libraries(DT3_benchmark_with_TBB PRIVATE benchmark::benchmark)

  add_executable(DT3_benchmark_with_TBB_CCC_approximate_size DT3_benchmark_with_TBB.cpp)
  CGAL_target_use_TBB(DT3_benchmark_with_TBB_CCC_approximate_size)
  target_compile_definitions(DT3_benchmark_with_TBB_CCC_approximate_size PRIVATE CGAL_CONCURRENT_COMPACT_CONTAINER_APPROXIMATE_SIZE)
  target_link_libraries(DT3_benchmark_with_TBB_CCC_approximate_size PRIVATE CGAL::CGAL benchmark::benchmark)
endif()
