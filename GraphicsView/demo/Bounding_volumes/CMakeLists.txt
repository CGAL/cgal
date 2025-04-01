# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Bounding_volumes_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets)

find_package(Eigen3 3.1.91 QUIET) #(requires 3.1.91 or greater)
include(CGAL_Eigen3_support)

if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This demo requires CGAL and Eigen, and will not be compiled.")
  return()
endif()

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  add_compile_definitions(QT_NO_KEYWORDS)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  qt_add_executable(
    Bounding_volumes Bounding_volumes.cpp Bounding_volumes.ui Bounding_volumes.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Bounding_volumes)

  target_link_libraries(Bounding_volumes PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                                 Qt6::Widgets
                                                 CGAL::Eigen3_support)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Bounding_volumes)

else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
