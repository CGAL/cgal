# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Segment_Delaunay_graph_Linf_2_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5 Core)

set(QT_USE_QTXML TRUE)
set(QT_USE_QTMAIN TRUE)
set(QT_USE_QTSCRIPT TRUE)
set(QT_USE_QTOPENGL TRUE)

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)
include_directories(BEFORE ./include)

if(CGAL_Qt5_FOUND AND Qt5_FOUND)
  set(CMAKE_AUTOMOC ON)
  include(${CGAL_USE_FILE})
  add_definitions(-DQT_NO_KEYWORDS)

  #--------------------------------
  # The "Segment Voronoi Linf" demo: Segment_voronoi_linf_2
  #--------------------------------
  # UI files (Qt Designer files)
  qt5_wrap_ui(CDT_UI_FILES Segment_voronoi_2.ui)

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources(CGAL_Qt5_RESOURCE_FILES ./Segment_voronoi_2.qrc)

  # The executable itself.
  add_executable(
    Segment_voronoi_linf_2 Segment_voronoi_linf_2.cpp ${CDT_UI_FILES}
                           ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Segment_voronoi_linf_2)

  target_link_libraries(Segment_voronoi_linf_2 PRIVATE CGAL::CGAL
                                                       CGAL::CGAL_Qt5 Qt5::Gui)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Segment_voronoi_linf_2)
else()

  message(
    STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
