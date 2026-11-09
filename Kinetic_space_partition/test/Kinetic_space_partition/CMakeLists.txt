# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)

project(Kinetic_space_partition_Tests)

find_package(CGAL REQUIRED)
include(CGAL_CreateSingleSourceCGALProgram)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

find_package(TBB QUIET)
include(CGAL_TBB_support)

set(targets kinetic_3d_test_all issue_8624)

set(project_linked_libraries)
set(project_compilation_definitions)

if(TARGET CGAL::TBB_support)
  message(STATUS "Found TBB")
  list(APPEND project_linked_libraries CGAL::TBB_support)
else()
  message(STATUS "NOTICE: Intel TBB was not found. Sequential code will be used.")
endif()

foreach(target ${targets})
  create_single_source_cgal_program("${target}.cpp")
  if(TARGET ${target})
    target_link_libraries(${target} PRIVATE ${project_linked_libraries} CGAL::Eigen3_support)
    target_compile_definitions(${target} PRIVATE ${project_compilation_definitions})
  endif()
endforeach()
