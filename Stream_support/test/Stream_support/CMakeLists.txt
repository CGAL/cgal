# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
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

if (NOT MSVC_VERSION OR MSVC_VERSION GREATER_EQUAL 1919 OR MSVC_VERSION LESS 1910)
  find_package(LASLIB QUIET)
  include(CGAL_LASLIB_support)
else()
  message(STATUS "NOTICE : the LAS reader does not work with your version of Visual Studio 2017.")
endif()


set(CMAKE_POLICY_DEFAULT_CMP0167 NEW)
find_package(VTK 9.0 QUIET COMPONENTS CommonCore IOCore IOLegacy IOXML FiltersCore FiltersSources)
if (VTK_FOUND AND VTK_LIBRARIES)
  message(STATUS "VTK ${VTK_VERSION} found ${VTK_LIBRARIES}")
else()
  message(STATUS "Tests that use VTK will not be compiled.")
endif() #VTK_FOUND

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
    if("${cppfile}" STREQUAL "test_LAS.cpp")
      if(TARGET CGAL::LASLIB_support)
        create_single_source_cgal_program("test_LAS.cpp")
        target_link_libraries(test_LAS PRIVATE CGAL::LASLIB_support)
      else()
        message(STATUS "NOTICE: Some tests require the LASlib library, and will not be compiled.")
      endif()
    else()
      if("${cppfile}" STREQUAL "test_VTK.cpp")
        if (VTK_FOUND AND VTK_LIBRARIES)
          create_single_source_cgal_program("test_VTK.cpp")
          target_link_libraries(test_VTK PRIVATE ${VTK_LIBRARIES})
          target_compile_definitions(test_VTK PRIVATE -DCGAL_USE_VTK -DNOMINMAX)
        endif()
      else()
        create_single_source_cgal_program("${cppfile}")
      endif()
    endif()
  endif()
endforeach()
