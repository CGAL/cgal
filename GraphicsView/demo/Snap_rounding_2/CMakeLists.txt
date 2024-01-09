# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Snap_rounding_2_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)
  add_definitions(-DQT_NO_KEYWORDS)

  qt_add_executable(
    Snap_rounding_2 Snap_rounding_2.cpp
    Snap_rounding_2.ui Snap_rounding_2.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Snap_rounding_2)

  target_link_libraries(Snap_rounding_2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                                Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Snap_rounding_2)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Snap_rounding_2)
else()
  message("NOTICE: This demo requires CGAL and Qt6<, and will not be compiled.")
endif()
