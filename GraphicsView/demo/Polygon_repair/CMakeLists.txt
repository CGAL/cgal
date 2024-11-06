# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Polygon_repair_Demo)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)
find_package(Qt6 QUIET COMPONENTS Widgets)

if(NOT CGAL_Qt6_FOUND OR NOT Qt6_FOUND)
  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")
  return()
endif()

add_definitions(-DQT_NO_KEYWORDS)
set(CMAKE_INCLUDE_CURRENT_DIR ON)


qt_add_executable(
  Boolean_2 Boolean_2.cpp
  Boolean_2.ui Boolean_2.qrc)
target_link_libraries(Boolean_2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                                      Qt6::Widgets)

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Boolean_2)

include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)

cgal_add_compilation_test(Boolean_2)

