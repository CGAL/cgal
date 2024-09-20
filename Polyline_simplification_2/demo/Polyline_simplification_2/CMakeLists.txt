# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Polyline_simplification_2_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt6 QUIET COMPONENTS Widgets)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  qt_add_executable(
    Polyline_simplification_2
    Polyline_simplification_2.cpp Polyline_simplification_2.ui Polyline_simplification_2.qrc)

  target_link_libraries(Polyline_simplification_2
                        PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Polyline_simplification_2)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Polyline_simplification_2)

else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
