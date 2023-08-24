# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Alpha_shapes_3_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets OpenGL)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  add_definitions(-DQT_NO_KEYWORDS)

  # ui file, created with Qt Designer
  qt6_wrap_ui(uis MainWindow.ui)

  # qrc files (resources files, that contain icons, at least)
  qt6_add_resources(CGAL_Qt6_RESOURCE_FILES ./Alpha_shape_3.qrc)

  add_executable(Alpha_shape_3 Alpha_shape_3.cpp MainWindow.cpp Viewer.cpp ${uis})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Alpha_shape_3)

  target_link_libraries(Alpha_shape_3 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                              Qt6::Widgets Qt6::OpenGL)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Alpha_shape_3)

else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
