cmake_minimum_required(VERSION 3.1...3.15)
project( benchmark_for_closest_rotation )



find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})

  find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
  if (EIGEN3_FOUND)
    create_single_source_cgal_program( "benchmark_for_concept_models.cpp" )
    CGAL_target_use_Eigen(benchmark_for_concept_models)
  else()
    message(STATUS "This program requires the Eigen library, version 3.1 or later and will not be compiled.")
  endif()
else()
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
endif()

