# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Largest_empty_rect_2_Demo)

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
  set(CMAKE_INCLUDE_CURRENT_DIR ON)
  set(CMAKE_AUTOMOC ON)

  #--------------------------------
  # Demo: Largest_empty_rectangle_2
  #--------------------------------
  # UI files (Qt Designer files)
  qt5_wrap_ui(DT_UI_FILES Largest_empty_rectangle_2.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./Largest_empty_rectangle_2.qrc)

  # The executable itself.
  add_executable(
    Largest_empty_rectangle_2 Largest_empty_rectangle_2.cpp ${DT_UI_FILES}
                              ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Largest_empty_rectangle_2)

  target_link_libraries(Largest_empty_rectangle_2
                        PRIVATE CGAL::CGAL CGAL::CGAL_Qt5 Qt5::Gui)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Largest_empty_rectangle_2)
else()

  message(
    STATUS "NOTICE: This demo requires CGAL, and Qt5, and will not be compiled."
  )

endif()
