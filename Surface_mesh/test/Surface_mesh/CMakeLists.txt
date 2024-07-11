# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Surface_mesh_Tests)

find_package(CGAL REQUIRED)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()

find_path(3MF_INCLUDE_DIR
  NAMES Model/COM/NMR_DLLInterfaces.h
  DOC "Path to lib3MF headers"
  )
find_library(3MF_LIBRARIES NAMES 3MF DOC "Path to the lib3MF library")
if(3MF_LIBRARIES AND 3MF_INCLUDE_DIR  AND EXISTS "${3MF_INCLUDE_DIR}/Model/COM/NMR_DLLInterfaces.h")
  include_directories(${3MF_INCLUDE_DIR})
  target_link_libraries(test_deprecated_io_sm PRIVATE ${3MF_LIBRARIES})
  target_compile_definitions(test_deprecated_io_sm PRIVATE -DCGAL_LINKED_WITH_3MF)
else()
  message(STATUS "NOTICE: read_3mf requires the lib3MF library, and will not be tested.")
endif()
