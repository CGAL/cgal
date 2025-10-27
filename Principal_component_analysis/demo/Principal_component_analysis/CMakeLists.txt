# This is the CMake script for compiling the PCA demo.

cmake_minimum_required(VERSION 3.12...3.29)
project(Principal_component_analysis_Demo)

include_directories(./)

# Find CGAL and CGAL Qt6
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(NOT TARGET CGAL::Eigen3_support)
  message("NOTICE: This project requires the Eigen library, and will not be compiled.")
  return()
endif()

# Find Qt6 itself
find_package(Qt6 QUIET COMPONENTS Widgets OpenGL)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  qt_add_executable(PCA_demo PCA_demo.cpp MainWindow.ui MainWindow.cpp Viewer.cpp Scene.cpp)

  target_link_libraries(PCA_demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                         CGAL::Eigen3_support Qt6::Widgets Qt6::OpenGL)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS PCA_demo)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(PCA_demo)

else(CGAL_Qt6_FOUND AND Qt6_FOUND)

  set(PCA_MISSING_DEPS "")

  if(NOT CGAL_Qt6_FOUND)
    set(PCA_MISSING_DEPS "the CGAL Qt6 library, ${PCA_MISSING_DEPS}")
  endif()

  if(NOT Qt6_FOUND)
    set(PCA_MISSING_DEPS "Qt6, ${PCA_MISSING_DEPS}")
  endif()

  message("NOTICE: This demo requires ${PCA_MISSING_DEPS} and will not be compiled.")

endif(CGAL_Qt6_FOUND AND Qt6_FOUND)
