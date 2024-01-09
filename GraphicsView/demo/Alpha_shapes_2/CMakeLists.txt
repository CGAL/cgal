# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Alpha_shapes_2_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  # The executable itself.
  qt_add_executable(
    Alpha_shapes_2 Alpha_shapes_2.cpp Alpha_shapes_2.ui Alpha_shapes_2.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Alpha_shapes_2)

  target_link_libraries(Alpha_shapes_2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                               Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Alpha_shapes_2)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Alpha_shapes_2)
else()
  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")
endif()
