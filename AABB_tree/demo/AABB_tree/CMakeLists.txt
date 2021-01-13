# This is the CMake script for compiling the AABB tree demo.

cmake_minimum_required(VERSION 3.1...3.15)
project(AABB_tree_Demo)

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

# Include this package's headers first
include_directories(BEFORE ./ ./include)

# Find CGAL and CGAL Qt5
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Gui Svg)

if(CGAL_Qt5_FOUND AND Qt5_FOUND)

  qt5_wrap_ui(UI_FILES MainWindow.ui)

  include(AddFileDependencies)

  qt5_generate_moc("MainWindow.h"
                   "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp")
  add_file_dependencies(MainWindow_moc.cpp
                        "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h")

  qt5_generate_moc("Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp")
  add_file_dependencies(Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h")

  qt5_generate_moc("Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp")
  add_file_dependencies(Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h")

  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES AABB_demo.qrc)

  add_file_dependencies(
    AABB_demo.cpp "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp")

  add_executable(
    AABB_demo AABB_demo.cpp ${UI_FILES} ${CGAL_Qt5_RESOURCE_FILES}
    #${CGAL_Qt5_MOC_FILES}
  )
  # Link with Qt libraries
  target_link_libraries(AABB_demo PRIVATE Qt5::OpenGL Qt5::Gui Qt5::Xml
                                          CGAL::CGAL CGAL::CGAL_Qt5)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS AABB_demo)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(AABB_demo)

else(CGAL_Qt5_FOUND
  AND Qt5_FOUND)

  set(AABB_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(AABB_MISSING_DEPS "CGAL_Qt5, ${AABB_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(AABB_MISSING_DEPS "Qt5, ${AABB_MISSING_DEPS}")
  endif()

  message(
    STATUS
      "NOTICE: This demo requires ${AABB_MISSING_DEPS}and will not be compiled."
  )

endif(
  CGAL_Qt5_FOUND
  AND Qt5_FOUND)
