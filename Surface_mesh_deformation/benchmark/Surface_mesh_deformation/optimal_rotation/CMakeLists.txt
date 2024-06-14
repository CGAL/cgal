cmake_minimum_required(VERSION 3.12...3.29)
project(benchmark_for_closest_rotation)

find_package(CGAL REQUIRED COMPONENTS Core)

find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("benchmark_for_concept_models.cpp")
  target_link_libraries(benchmark_for_concept_models PUBLIC CGAL::Eigen3_support)
else()
  message("NOTICE: This program requires requires Eigen 3.1 (or greater) or later and will not be compiled.")
endif()
