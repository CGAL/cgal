# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Segment_Delaunay_graph_Linf_2_Demo)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6 Core)

find_package(Qt6 QUIET COMPONENTS Widgets)
include_directories(BEFORE ./include)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)
  add_compile_definitions(QT_NO_KEYWORDS)

  qt_add_executable(
    Segment_voronoi_linf_2 Segment_voronoi_linf_2.cpp
    Segment_voronoi_2.ui Segment_voronoi_2.qrc)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Segment_voronoi_linf_2)

  target_link_libraries(Segment_voronoi_linf_2 PRIVATE CGAL::CGAL
                                                       CGAL::CGAL_Qt6 Qt6::Widgets)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Segment_voronoi_linf_2)
else()
  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")
endif()
