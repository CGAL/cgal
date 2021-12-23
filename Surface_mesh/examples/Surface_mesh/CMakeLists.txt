# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

#/!\ /!\
#/!\ /!\
# Used in /CGAL/Documentation/doc/Documentation/Developer_manual/create_and_use_a_cmakelist.txt.
# Careful when modifying

cmake_minimum_required(VERSION 3.1...3.22)
project(Surface_mesh_Examples)

#CGAL_Qt5 is needed for the drawing.
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

create_single_source_cgal_program("sm_derivation.cpp")
create_single_source_cgal_program("sm_join.cpp")
create_single_source_cgal_program("sm_aabbtree.cpp")
create_single_source_cgal_program("sm_bgl.cpp")
create_single_source_cgal_program("sm_circulators.cpp")
create_single_source_cgal_program("sm_do_intersect.cpp")
create_single_source_cgal_program("sm_iterators.cpp")
create_single_source_cgal_program("sm_kruskal.cpp")
create_single_source_cgal_program("sm_memory.cpp")
create_single_source_cgal_program("sm_properties.cpp")
create_single_source_cgal_program("sm_draw_small_faces.cpp")
create_single_source_cgal_program("check_orientation.cpp")

#create the executable of the application
create_single_source_cgal_program("draw_surface_mesh.cpp")
if(CGAL_Qt5_FOUND)

  #link it with the required CGAL libraries
  target_link_libraries(draw_surface_mesh PUBLIC CGAL::CGAL_Basic_viewer)

  target_link_libraries(sm_draw_small_faces PUBLIC CGAL::CGAL_Basic_viewer)
endif()

#end of the file
