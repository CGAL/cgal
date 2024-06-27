# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Generator_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("random_hull_test.cpp")
create_single_source_cgal_program("random_poly_test.cpp")
create_single_source_cgal_program("rcs_test.cpp")
create_single_source_cgal_program("test_combination_enumerator.cpp")
create_single_source_cgal_program("test_generators.cpp")
create_single_source_cgal_program("test_tetrahedron_3.cpp")
create_single_source_cgal_program("test_triangle_2.cpp")
create_single_source_cgal_program("test_triangle_3.cpp")

find_package(Eigen3 3.1.0 QUIET) #(3.1.0 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("generic_random_test.cpp")
  target_link_libraries(generic_random_test PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: The test 'generic_random_test' uses Eigen, and will not be compiled.")
endif()
