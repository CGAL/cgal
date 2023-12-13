# Created by the script cgal_create_CMakeLists.
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.15)

project(Kinetic_space_partition_Examples)

set(CMAKE_CXX_STANDARD 14)

find_package(CGAL REQUIRED COMPONENTS Core)
include(CGAL_CreateSingleSourceCGALProgram)

find_package(Boost REQUIRED)
if(Boost_FOUND)
  message(STATUS "Found Boost")

  find_package(Eigen3 3.1.0 REQUIRED)
  if(Eigen3_FOUND)
    message(STATUS "Found Eigen")
    include(CGAL_Eigen_support)

    set(targets kinetic_partition)

    set(project_linked_libraries)
    set(project_compilation_definitions)

    foreach(target ${targets})
      create_single_source_cgal_program("${target}.cpp")
      if(TARGET ${target})
        target_link_libraries(${target} PUBLIC ${project_linked_libraries} CGAL::Eigen_support)
        target_compile_definitions(${target} PUBLIC ${project_compilation_definitions})
      endif()
    endforeach()
  else()
    message(ERROR "This program requires the Eigen library, and will not be compiled.")
  endif()
else()
  message(ERROR "This program requires the Boost library, and will not be compiled.")
endif()
