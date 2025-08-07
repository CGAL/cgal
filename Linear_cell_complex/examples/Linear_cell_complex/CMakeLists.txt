# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Linear_cell_complex_Examples)

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# For Gprof.
# ADD_DEFINITIONS("-pg")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# To use valgrind, we must disable rounding math check.
# add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK)
# create_single_source_cgal_program("hexmeshing/hexmeshing_unit_test.cpp")
# create_single_source_cgal_program("hexmeshing/hexmeshing_simple_test.cpp")
create_single_source_cgal_program("hexmeshing/hexmeshing_demo.cpp")
# create_single_source_cgal_program("hexmeshing/hexmeshing_tests.cpp")

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

create_single_source_cgal_program("query_replace/hexa-subdivision.cpp")
create_single_source_cgal_program("query_replace/tetra-to-hexa.cpp")


if(CGAL_Qt6_FOUND)
  target_link_libraries(draw_linear_cell_complex PUBLIC CGAL::CGAL_Basic_viewer)
  # target_link_libraries(hexmeshing_unit_test PUBLIC CGAL::CGAL_Basic_viewer)
  # target_link_libraries(hexmeshing_simple_test PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(hexmeshing_demo PUBLIC CGAL::CGAL_Basic_viewer)
  # target_link_libraries(hexmeshing_tests PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(linear_cell_complex_3_incremental_builder PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(linear_cell_complex_3_insert PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(hexa-subdivision PUBLIC CGAL::CGAL_Basic_viewer)
  target_link_libraries(tetra-to-hexa PUBLIC CGAL::CGAL_Basic_viewer)
endif()

# target_link_libraries(hexmeshing_unit_test PUBLIC Eigen3::Eigen)
# target_link_libraries(hexmeshing_simple_test PUBLIC Eigen3::Eigen)
target_link_libraries(hexmeshing_demo PUBLIC Eigen3::Eigen)

# For now we use TetGen library, but later we'll remove the dependencies
find_library(TETGEN_LIBRARY libtet.a PATHS "/usr/local/lib")
find_path(TETGEN_INCLUDE_DIR tetgen.h PATHS "/usr/local/include")

if(NOT TETGEN_LIBRARY OR NOT TETGEN_INCLUDE_DIR)
  message(FATAL_ERROR "TetGen library or header not found")
endif()

target_include_directories(tetra-to-hexa PRIVATE ${TETGEN_INCLUDE_DIR})
target_link_libraries(tetra-to-hexa PRIVATE ${TETGEN_LIBRARY})