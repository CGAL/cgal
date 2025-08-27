# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Arrangement_on_surface_2_Examples)

find_package(CGAL REQUIRED COMPONENTS Core OPTIONAL_COMPONENTS Qt6)

# create a target per cppfile
file(
  GLOB cppfiles
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
foreach(cppfile ${cppfiles})
  create_single_source_cgal_program("${cppfile}")
endforeach()

if(CGAL_Qt6_FOUND)
  target_link_libraries(draw_arr PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(linear_conics PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(parabolas PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(ellipses PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(hyperbolas PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(polylines PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(circles PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(circular_arcs PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(spherical_insert PRIVATE CGAL::CGAL_Basic_viewer)
else()
  message(
    STATUS
      "NOTICE: Several examples require Qt and will not be compiled."
  )
endif()
