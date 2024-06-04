# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.14)
project( Polyline_distance_Examples )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )


find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)

create_single_source_cgal_program( "Frechet_distance_2.cpp" )

if(TARGET CGAL::Eigen3_support)
create_single_source_cgal_program( "Frechet_distance_d.cpp" )
target_link_libraries(Frechet_distance_d PUBLIC CGAL::Eigen3_support)
endif()