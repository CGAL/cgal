cmake_minimum_required(VERSION 3.12...3.31)
project( classical_Frechet_distance )

find_package(CGAL REQUIRED)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

create_single_source_cgal_program("Compute_classical_Frechet_distance_3.cpp")
create_single_source_cgal_program("txt2wkt_3.cpp")

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program( "Compute_classical_Frechet_distance_100.cpp" )
  target_link_libraries(Compute_classical_Frechet_distance_100 PUBLIC CGAL::Eigen3_support)
endif()
