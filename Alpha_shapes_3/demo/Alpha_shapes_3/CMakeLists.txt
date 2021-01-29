# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Alpha_shapes_3_Demo)

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

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

if(CGAL_Qt5_FOUND AND Qt5_FOUND)

  add_definitions(-DQT_NO_KEYWORDS)

  # include(${QT_USE_FILE})
  include_directories(BEFORE ./)

  # ui file, created wih Qt Designer
  qt5_wrap_ui(uis MainWindow.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./Alpha_shape_3.qrc)

  add_executable(
    Alpha_shape_3 Alpha_shape_3.cpp MainWindow.cpp Viewer.cpp ${uis}
                  ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Alpha_shape_3)

  target_link_libraries(Alpha_shape_3 PRIVATE CGAL::CGAL CGAL::CGAL_Qt5
                                              Qt5::OpenGL Qt5::Gui)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Alpha_shape_3)

else()

  message(
    STATUS "NOTICE: This demo requires CGAL, and Qt5, and will not be compiled."
  )

endif()
