# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Bounding_volumes_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets)
if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  qt_add_executable(
    Bounding_volumes Bounding_volumes.cpp Bounding_volumes.ui Bounding_volumes.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Bounding_volumes)

  target_link_libraries(Bounding_volumes PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                                 Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Bounding_volumes)

else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
