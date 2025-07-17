# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project( Frechet_distance_Examples )

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core )


find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)

create_single_source_cgal_program( "Frechet_distance_2.cpp" )
create_single_source_cgal_program( "Frechet_distance_3.cpp" )
create_single_source_cgal_program( "Frechet_DS_2.cpp" )
create_single_source_cgal_program( "Frechet_DS_3.cpp" )

if(TARGET CGAL::Eigen3_support)
create_single_source_cgal_program( "Frechet_DS_d.cpp" )
target_link_libraries(Frechet_DS_d PUBLIC CGAL::Eigen3_support)
create_single_source_cgal_program( "Frechet_distance_d_2.cpp" )
target_link_libraries(Frechet_distance_d_2 PUBLIC CGAL::Eigen3_support)

create_single_source_cgal_program( "Frechet_distance_d.cpp" )
target_link_libraries(Frechet_distance_d PUBLIC CGAL::Eigen3_support)
endif()


create_single_source_cgal_program("parallel_Frechet_DS_2.cpp")
find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  message(STATUS "Found TBB")

  target_link_libraries(parallel_Frechet_DS_2 PRIVATE CGAL::TBB_support)
else()
  message(STATUS "NOTICE: The example 'parallel_Frechet_DS_2' requires TBB, and will not be compiled")
endif()
