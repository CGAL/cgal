# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)

project (Triangulation_on_sphere_2_Demo)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()


# Find CGAL and CGAL Qt5
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Script OpenGL Gui Svg)

find_package(Eigen3 3.1.0) #(requires 3.1.0 or greater)
include(CGAL_Eigen3_support)

if(CGAL_Qt5_FOUND AND Qt5_FOUND AND TARGET CGAL::Eigen3_support)
  # Include this package's headers first

  include_directories(BEFORE ./ ./include)
  # ui file, created wih Qt Designer
  qt5_wrap_ui( uis Mainwindow.ui )

  #qt5_generate_moc( main.cpp Mainwindow.moc)

  add_executable  ( Triangulation_on_sphere_2_Demo main.cpp Viewer.cpp ${uis})

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Triangulation_on_sphere_2_Demo )

  target_link_libraries( Triangulation_on_sphere_2_Demo PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 CGAL::Eigen3_support)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test( Triangulation_on_sphere_2_Demo )
else()
  message(STATUS "NOTICE: This demo requires CGAL, Qt5 and OpenGL, and will not be compiled.")
endif()
