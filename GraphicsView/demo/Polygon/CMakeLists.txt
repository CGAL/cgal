# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Polygon_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6 Core)

find_package(Eigen3 3.1.0 QUIET) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)
if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This demo requires the Eigen library, and will not be compiled.")
  return()
endif()

find_package(Qt6 QUIET COMPONENTS Widgets)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  if(CGAL_Core_FOUND)
    add_definitions(-DCGAL_USE_CORE)
  endif()

  qt_add_executable(Polygon_2 Polygon_2.cpp
                    Polygon_2.ui Polygon_2.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Polygon_2)

  target_link_libraries(Polygon_2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                          CGAL::Eigen3_support Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Polygon_2)
else()
  message("NOTICE: This demo requires CGAL, CGAL_Core, and Qt6, and will not be compiled.")
endif()
