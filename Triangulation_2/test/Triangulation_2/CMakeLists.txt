# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Triangulation_2_Tests)

find_package(CGAL REQUIRED)

include_directories(BEFORE "include")

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()

if(CGAL_ENABLE_TESTING)
  set_tests_properties(
    "execution   of  test_constrained_triangulation_2"
    "execution   of  test_delaunay_triangulation_2"
    "execution   of  test_triangulation_geom_traits"
    "execution   of  test_triangulation_2"
    "execution   of  test_triangulation_2_bis"
    "execution   of  test_delaunay_hierarchy_2"
    "execution   of  test_const_triang_plus_2"
    "execution   of  test_regular_triangulation_2"
    "execution   of  test_const_del_triangulation_2"
    "execution   of  test_regular_hierarchy_2"
    "execution   of  test_deprecated_projection_traits"
    PROPERTIES RESOURCE_LOCK Triangulation_2_Tests_IO)
endif()
