# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.12...3.31)
project(Generalized_map_Tests)

# CGAL and its components
find_package(CGAL REQUIRED)

set(hfiles Generalized_map_2_test.h
           Generalized_map_3_test.h
           Generalized_map_4_test.h
           GMap_test_insertions.h)

create_single_source_cgal_program("Generalized_map_test.cpp" ${hfiles})

add_executable(Generalized_map_test_index  Generalized_map_test.cpp ${hfiles})
target_compile_definitions(Generalized_map_test_index PRIVATE USE_COMPACT_CONTAINER_WITH_INDEX)
target_link_libraries(Generalized_map_test_index PRIVATE CGAL CGAL::Data)
cgal_add_compilation_test(Generalized_map_test_index)

create_single_source_cgal_program("gmap_test_split_attribute.cpp")
