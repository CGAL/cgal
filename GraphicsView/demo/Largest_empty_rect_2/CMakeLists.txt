# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Largest_empty_rect_2_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  add_compile_definitions(QT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)


  qt_add_executable(
    Largest_empty_rectangle_2 Largest_empty_rectangle_2.cpp
                              Largest_empty_rectangle_2.ui Largest_empty_rectangle_2.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Largest_empty_rectangle_2)

  target_link_libraries(Largest_empty_rectangle_2
                        PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Largest_empty_rectangle_2)
else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
