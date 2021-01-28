# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.15)
project(Linear_cell_complex_Examples)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

if(CGAL_Qt5_FOUND)
  add_definitions(-DCGAL_USE_BASIC_VIEWER -DQT_NO_KEYWORDS)
endif()

# For Gprof.
# ADD_DEFINITIONS("-pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# To use valgrind, we must disable rounding math ckeck.
# add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)

create_single_source_cgal_program("gmap_linear_cell_complex_3.cpp")
create_single_source_cgal_program("linear_cell_complex_3.cpp")
create_single_source_cgal_program(
  "linear_cell_complex_3_attributes_management.cpp")
create_single_source_cgal_program("linear_cell_complex_3_operations.cpp")
create_single_source_cgal_program(
  "linear_cell_complex_3_with_colored_vertices.cpp")
create_single_source_cgal_program("linear_cell_complex_3_with_mypoint.cpp")
create_single_source_cgal_program("linear_cell_complex_4.cpp")
create_single_source_cgal_program("plane_graph_to_lcc_2.cpp")
create_single_source_cgal_program("voronoi_2.cpp")
create_single_source_cgal_program("voronoi_3.cpp")

create_single_source_cgal_program("draw_linear_cell_complex.cpp")
if(CGAL_Qt5_FOUND)
  target_link_libraries(draw_linear_cell_complex PUBLIC CGAL::CGAL_Qt5)
endif()
