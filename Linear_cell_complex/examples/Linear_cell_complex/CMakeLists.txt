# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Linear_cell_complex_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)
find_package(Eigen3 3.3 QUIET NO_MODULE)

# For Gprof.
# ADD_DEFINITIONS("-pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# To use valgrind, we must disable rounding math check.
# add_compile_definitions(CGAL_DISABLE_ROUNDING_MATH_CHECK)

create_single_source_cgal_program("gmap_linear_cell_complex_3.cpp")
create_single_source_cgal_program("linear_cell_complex_3.cpp")
create_single_source_cgal_program(
  "linear_cell_complex_3_attributes_management.cpp")
create_single_source_cgal_program("linear_cell_complex_3_incremental_builder.cpp")
create_single_source_cgal_program("linear_cell_complex_3_insert.cpp")
create_single_source_cgal_program("linear_cell_complex_3_operations.cpp")
create_single_source_cgal_program(
  "linear_cell_complex_3_with_colored_vertices.cpp")
create_single_source_cgal_program("linear_cell_complex_3_with_mypoint.cpp")
create_single_source_cgal_program("linear_cell_complex_4.cpp")
create_single_source_cgal_program("plane_graph_to_lcc_2.cpp")
create_single_source_cgal_program("voronoi_2.cpp")
create_single_source_cgal_program("voronoi_3.cpp")

create_single_source_cgal_program("draw_linear_cell_complex.cpp")

if(CGAL_Qt6_FOUND)
  target_link_libraries(draw_linear_cell_complex PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(linear_cell_complex_3_incremental_builder PRIVATE CGAL::CGAL_Basic_viewer)
  target_link_libraries(linear_cell_complex_3_insert PRIVATE CGAL::CGAL_Basic_viewer)
endif()

if(Eigen3_FOUND)
  message(STATUS "Eigen3 found: enabling hexmeshing_* targets")

  create_single_source_cgal_program("hexmeshing_cube.cpp")

  target_link_libraries(hexmeshing_cube PRIVATE Eigen3::Eigen)

  if(CGAL_Qt6_FOUND)
    target_link_libraries(hexmeshing_cube PRIVATE CGAL::CGAL_Basic_viewer)
  endif()
else()
  message(STATUS "Eigen3 NOT found: hexmeshing_* targets will be skipped")
endif()
