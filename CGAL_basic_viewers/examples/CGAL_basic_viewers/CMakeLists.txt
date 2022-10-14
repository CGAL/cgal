# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

#/!\ /!\
#/!\ /!\
# Used in /CGAL/Documentation/doc/Documentation/Developer_manual/create_and_use_a_cmakelist.txt.
# Careful when modifying

cmake_minimum_required(VERSION 3.1...3.23)
project(CGAL_basic_viewers)

#CGAL_Qt5 is needed for the drawing.
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

create_single_source_cgal_program("draw_surface_mesh_small_faces.cpp")
if(CGAL_Qt5_FOUND)
  #link it with the required CGAL libraries
  target_link_libraries(draw_surface_mesh_small_faces PUBLIC CGAL::CGAL_Basic_viewer)
endif()

#end of the file
