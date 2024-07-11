# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

#/!\ /!\
#/!\ /!\
# Used in /CGAL/Documentation/doc/Documentation/Developer_manual/create_and_use_a_cmakelist.txt.
# Careful when modifying

cmake_minimum_required(VERSION 3.12...3.29)
project(Basic_viewer_Examples)

#CGAL_Qt6 is needed for the drawing.
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)
find_package(Eigen3 3.1.0)
include(CGAL_Eigen3_support)

create_single_source_cgal_program("draw_lcc.cpp")
create_single_source_cgal_program("draw_mesh_and_points.cpp")
create_single_source_cgal_program("draw_several_windows.cpp")
create_single_source_cgal_program("draw_surface_mesh_height.cpp")
create_single_source_cgal_program("draw_surface_mesh_small_faces.cpp")
create_single_source_cgal_program("draw_surface_mesh_vcolor.cpp")

if(CGAL_Qt6_FOUND)
  #link it with the required CGAL libraries
  target_link_libraries(draw_lcc PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(draw_mesh_and_points PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(draw_several_windows PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(draw_surface_mesh_height PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(draw_surface_mesh_small_faces PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(draw_surface_mesh_vcolor PUBLIC CGAL::CGAL_Basic_viewer)
else()
  message("CGAL_Qt6 not configured: examples that require Qt will not be compiled.")
endif()

if(TARGET CGAL::Eigen3_support)
  target_link_libraries(draw_mesh_and_points PUBLIC CGAL::Eigen3_support)
  target_link_libraries(draw_several_windows PUBLIC CGAL::Eigen3_support)
endif()

#end of the file
