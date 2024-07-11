# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Nef_3_Tests)

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
    "execution   of  Nef_3_problematic_construction"
    "execution   of  test_nef_3_io_Cartesian"
    "execution   of  test_nef_3_io_Cartesian_Lazy"
    "execution   of  test_nef_3_io_EPEC"
    "execution   of  test_nef_3_io_EPIC"
    "execution   of  test_nef_3_io_Homogeneous"
    "execution   of  test_nef_3_io_Homogenoeus_Lazy"
    "execution   of  test_with_extended_homogeneous"
    "execution   of  test_with_homogeneous"
    "execution   of  test_with_lazy"
    PROPERTIES RESOURCE_LOCK Nef_3_Tests_IO)
endif()
