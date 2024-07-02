# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Alpha_shapes_3_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets OpenGL)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  add_definitions(-DQT_NO_KEYWORDS)

  # Instruct CMake to run moc/ui/rcc automatically when needed.
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  qt_add_executable(Alpha_shape_3 Alpha_shape_3.cpp MainWindow.cpp Viewer.cpp
                    MainWindow.ui Alpha_shape_3.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Alpha_shape_3)

  target_link_libraries(Alpha_shape_3 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                              Qt6::Widgets Qt6::OpenGL)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Alpha_shape_3)

else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
