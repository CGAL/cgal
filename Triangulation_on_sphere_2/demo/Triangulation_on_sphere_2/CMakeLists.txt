# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)

project (Triangulation_on_sphere_2_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# Find CGAL and CGAL Qt6
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

# Find Qt6 itself
find_package(Qt6 QUIET COMPONENTS Widgets OpenGL)

find_package(Eigen3 QUIET)
include(CGAL_Eigen3_support)

if(CGAL_Qt6_FOUND AND Qt6_FOUND AND TARGET CGAL::Eigen3_support)
  # Include this package's headers first
  include_directories(BEFORE ./ ./include)

  set(CMAKE_AUTOUIC ON)

  qt_add_executable  ( Triangulation_on_sphere_2_Demo main.cpp Viewer.cpp)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Triangulation_on_sphere_2_Demo )

  target_link_libraries( Triangulation_on_sphere_2_Demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 Qt6::OpenGL Qt6::Widgets CGAL::Eigen3_support)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test( Triangulation_on_sphere_2_Demo )
else()
  message("NOTICE: This demo requires CGAL, Qt6, and Eigen, and will not be compiled.")
endif()
