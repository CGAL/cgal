# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Snap_rounding_2_Demo)

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
  set(CMAKE_AUTOMOC ON)
  add_definitions(-DQT_NO_KEYWORDS)
  set(CMAKE_INCLUDE_CURRENT_DIR ON)

  #--------------------------------
  # UI files (Qt Designer files)
  qt5_wrap_ui(DT_UI_FILES Snap_rounding_2.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./Snap_rounding_2.qrc)

  # The executable itself.
  add_executable(
    Snap_rounding_2 Snap_rounding_2.cpp ${DT_UI_FILES}
                    ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Snap_rounding_2)

  target_link_libraries(Snap_rounding_2 PRIVATE CGAL::CGAL CGAL::CGAL_Qt5
                                                Qt5::Gui)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Snap_rounding_2)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Snap_rounding_2)

else()

  message(
    STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
