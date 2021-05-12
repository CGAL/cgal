cmake_minimum_required(VERSION 3.1...3.20)
project(benchmark_for_closest_rotation)

find_package(CGAL REQUIRED COMPONENTS Core)

include(${CGAL_USE_FILE})

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("benchmark_for_concept_models.cpp")
  target_link_libraries(benchmark_for_concept_models PUBLIC CGAL::Eigen3_support)
else()
  message(
    STATUS
      "This program requires the Eigen library, version 3.1 or later and will not be compiled."
  )
endif()
