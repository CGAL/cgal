# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Polygon_repair_Tests)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

# create a target per cppfile
#file(
#  GLOB cppfiles
#  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
#  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
#foreach(cppfile ${cppfiles})
#  create_single_source_cgal_program("${cppfile}")
#endforeach()

create_single_source_cgal_program( "draw_test_polygons.cpp" )
create_single_source_cgal_program( "exact_test.cpp")
create_single_source_cgal_program( "repair_polygon_2_test.cpp" )

if(CGAL_Qt6_FOUND)
  target_link_libraries(draw_test_polygons PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(exact_test PUBLIC CGAL::CGAL_Basic_viewer)
endif()
