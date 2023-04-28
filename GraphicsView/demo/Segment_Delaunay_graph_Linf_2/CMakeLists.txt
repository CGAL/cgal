# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Segment_Delaunay_graph_Linf_2_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6 Core)

set(QT_USE_QTXML TRUE)
set(QT_USE_QTMAIN TRUE)
set(QT_USE_QTSCRIPT TRUE)
set(QT_USE_QTOPENGL TRUE)

find_package(Qt6 QUIET COMPONENTS Widgets OpenGLWidgets)
include_directories(BEFORE ./include)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  set(CMAKE_AUTOMOC ON)
  include(${CGAL_USE_FILE})
  add_definitions(-DQT_NO_KEYWORDS)

  # UI files (Qt Designer files)
  qt6_wrap_ui(CDT_UI_FILES Segment_voronoi_2.ui)

  # qrc files (resources files, that contain icons, at least)
  qt6_add_resources(CGAL_Qt6_RESOURCE_FILES ./Segment_voronoi_2.qrc)

  # The executable itself.
  add_executable(
    Segment_voronoi_linf_2 Segment_voronoi_linf_2.cpp ${CDT_UI_FILES}
                           ${CGAL_Qt6_RESOURCE_FILES} ${CGAL_Qt6_MOC_FILES})

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Segment_voronoi_linf_2)

  target_link_libraries(Segment_voronoi_linf_2 PRIVATE CGAL::CGAL
                                                       CGAL::CGAL_Qt6 Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Segment_voronoi_linf_2)
else()
  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")
endif()
