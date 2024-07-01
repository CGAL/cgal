# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project( Isosurfacing_3_benchmark )

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)

find_package(TBB)
include(CGAL_TBB_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("benchmark.cpp" )
  create_single_source_cgal_program("contouring_seq_vs_parallel_implicit.cpp" )
  create_single_source_cgal_program("contouring_seq_vs_parallel_image.cpp" )

  target_compile_definitions(benchmark PUBLIC ${SCENARIO})
  target_compile_definitions(benchmark PUBLIC ${KERNEL})
  target_compile_definitions(benchmark PUBLIC ${ALGO})
  target_compile_definitions(benchmark PUBLIC ${TAG})

  target_link_libraries(benchmark PRIVATE CGAL::Eigen3_support)
  target_link_libraries(contouring_seq_vs_parallel_implicit PRIVATE CGAL::Eigen3_support)
  target_link_libraries(contouring_seq_vs_parallel_image PRIVATE CGAL::Eigen3_support)

  if(TARGET CGAL::TBB_support)
    target_link_libraries(benchmark PUBLIC CGAL::TBB_support)
    target_link_libraries(contouring_seq_vs_parallel_implicit PUBLIC CGAL::TBB_support)
    target_link_libraries(contouring_seq_vs_parallel_image PUBLIC CGAL::TBB_support)
  endif()
else()
  message(STATUS "NOTICE: Some benchmarks use Eigen, and will not be compiled.")
endif()
