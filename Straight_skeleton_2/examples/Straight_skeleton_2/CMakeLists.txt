# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project( Straight_skeleton_2_Examples )

find_package(CGAL REQUIRED COMPONENTS Qt6 Core)

include_directories(BEFORE include)

# create a target per cppfile
file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program( "${cppfile}" )
endforeach()

if(CGAL_Qt6_FOUND)
  target_link_libraries(draw_straight_skeleton_2 PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(exterior_offset_of_multiple_polygons_with_holes PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(Create_straight_skeleton_2 PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(Show_offset_polygon PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(Show_straight_skeleton PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(extrude_skeleton PRIVATE CGAL::CGAL_Basic_viewer)
else()
  message(STATUS "NOTICE: The example 'draw_straight_skeleton_2' requires Qt and will not be compiled.")
endif()
