# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Triangulation_2_Demo)

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

#--------------------------------
# The "constrained Delaunay" demo: Constrained_Delaunay_triangulation_2
#--------------------------------

qt_add_executable(
  Constrained_Delaunay_triangulation_2
  Constrained_Delaunay_triangulation_2.cpp
  Constrained_Delaunay_triangulation_2.ui Constrained_Delaunay_triangulation_2.qrc)
target_link_libraries(Constrained_Delaunay_triangulation_2
                      PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)
target_include_directories(Constrained_Delaunay_triangulation_2
                           PRIVATE ./include)

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Constrained_Delaunay_triangulation_2)

#--------------------------------
# The "Delaunay" demo: Delaunay_triangulation_2
#--------------------------------
qt_add_executable(
  Delaunay_triangulation_2 Delaunay_triangulation_2.cpp
  Delaunay_triangulation_2.ui Delaunay_triangulation_2.qrc)
target_link_libraries(Delaunay_triangulation_2
                      PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Delaunay_triangulation_2)

#--------------------------------
# The "Regular" demo: Regular_triangulation_2
#--------------------------------

qt_add_executable(
  Regular_triangulation_2 Regular_triangulation_2.cpp
  Regular_triangulation_2.ui Regular_triangulation_2.qrc)
target_link_libraries(Regular_triangulation_2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                                      Qt6::Widgets)

add_to_cached_list(CGAL_EXECUTABLE_TARGETS Regular_triangulation_2)

include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
foreach(target Constrained_Delaunay_triangulation_2 Delaunay_triangulation_2
               Regular_triangulation_2)
  cgal_add_compilation_test(${target})
endforeach()
