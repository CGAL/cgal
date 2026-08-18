# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Vector_graphics_on_surfaces_Tests)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)

include(CGAL_Eigen3_support)

create_single_source_cgal_program( "test_straightest_geodesic.cpp" )
target_link_libraries(test_straightest_geodesic PRIVATE CGAL::Eigen3_support)
create_single_source_cgal_program("test_Bsurf.cpp")
target_link_libraries(test_Bsurf PUBLIC CGAL::Eigen3_support)
create_single_source_cgal_program("test_Bsurf_grid.cpp")
target_link_libraries(test_Bsurf_grid PUBLIC CGAL::Eigen3_support)
