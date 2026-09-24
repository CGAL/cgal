# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Spatial_searching_Tests)

find_package(CGAL REQUIRED)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

list(REMOVE_ITEM cppfiles "issue_9534.cpp")

foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("issue_9534.cpp")
  target_link_libraries(issue_9534 PRIVATE CGAL::Eigen3_support)
else()
  message(STATUS "NOTICE: The test 'issue_9534' requires Eigen3 and will not be compiled.")
endif()

find_package(TBB QUIET)
include(CGAL_TBB_support)
if(TARGET CGAL::TBB_support)
  target_link_libraries(test_parallel_build PRIVATE CGAL::TBB_support)
  if(CGAL_ENABLE_TESTING)
    set_property(TEST "execution   of  test_parallel_build"
      PROPERTY RUN_SERIAL 1)
  endif()
else()
  message(STATUS "NOTICE: Intel TBB was not found. Sequential code will be used.")
endif()
