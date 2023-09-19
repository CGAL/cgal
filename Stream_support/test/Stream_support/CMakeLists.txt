# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Stream_support_Tests)
find_package(CGAL REQUIRED)
find_path(
  3MF_INCLUDE_DIR
  NAMES Model/COM/NMR_DLLInterfaces.h
  DOC "Path to lib3MF headers")
find_library(
  3MF_LIBRARIES
  NAMES 3MF
  DOC "Path to the lib3MF library")

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  if("${cppfile}" STREQUAL "test_3mf_to_sm.cpp")
    if(3MF_LIBRARIES
       AND 3MF_INCLUDE_DIR
       AND EXISTS "${3MF_INCLUDE_DIR}/Model/COM/NMR_DLLInterfaces.h")
      include_directories(${3MF_INCLUDE_DIR})
      create_single_source_cgal_program("${cppfile}")
      target_link_libraries(test_3mf_to_sm PRIVATE ${3MF_LIBRARIES})
    else()
      message(STATUS "NOTICE: Some tests require the lib3MF library, and will not be compiled.")
    endif()
  else()
    create_single_source_cgal_program("${cppfile}")
  endif()
endforeach()
