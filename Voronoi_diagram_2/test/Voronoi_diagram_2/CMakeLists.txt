# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.22)
project(Voronoi_diagram_2_Tests)

find_package(CGAL REQUIRED)

include_directories(BEFORE "include")

create_single_source_cgal_program("vda_ag.cpp")
create_single_source_cgal_program("vda_dt.cpp")
create_single_source_cgal_program("vda_pt.cpp")
create_single_source_cgal_program("vda_rt.cpp")
create_single_source_cgal_program("vda_sdg.cpp")

find_package(Eigen3 3.2.0) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("vda_tos2.cpp")
  target_link_libraries(vda_tos2 PUBLIC CGAL::Eigen3_support)
endif()
