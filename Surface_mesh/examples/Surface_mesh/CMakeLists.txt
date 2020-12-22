# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Surface_mesh_Examples)

if(POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  if(NOT POLICY CMP0070)
    cmake_policy(SET CMP0053 OLD)
  else()
    cmake_policy(SET CMP0053 NEW)
  endif()
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

#CGAL_Qt5 is needed for the drawing.
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

if(CGAL_Qt5_FOUND)
  #required to use basic_viewer
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

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

#create the executable of the application

create_single_source_cgal_program("draw_surface_mesh.cpp")
create_single_source_cgal_program("sm_draw_small_faces.cpp")
create_single_source_cgal_program("check_orientation.cpp")
if(CGAL_Qt5_FOUND)

  #link it with the required CGAL libraries

  target_link_libraries(draw_surface_mesh PUBLIC CGAL::CGAL_Qt5)
  target_link_libraries(sm_draw_small_faces PUBLIC CGAL::CGAL_Qt5)
endif()

#end of the file
