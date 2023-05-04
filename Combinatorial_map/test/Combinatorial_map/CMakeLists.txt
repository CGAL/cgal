# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Combinatorial_map_Tests)

find_package(CGAL REQUIRED)

set(hfiles Combinatorial_map_2_test.h Combinatorial_map_3_test.h
  Combinatorial_map_test_iterators.h Face_graph_wrapper_test.h)

# create a target per cppfile
create_single_source_cgal_program(Combinatorial_map_test.cpp ${hfiles})
create_single_source_cgal_program(Combinatorial_map_copy_test.cpp ${hfiles})

# Same targets, defining USE_COMPACT_CONTAINER_WITH_INDEX to test index version
add_executable(Combinatorial_map_test_index  Combinatorial_map_test.cpp ${hfiles})
target_compile_definitions(Combinatorial_map_test_index PUBLIC USE_COMPACT_CONTAINER_WITH_INDEX)
target_link_libraries(Combinatorial_map_test_index PUBLIC CGAL CGAL::Data)
cgal_add_compilation_test(Combinatorial_map_test_index)

add_executable(Combinatorial_map_copy_test_index Combinatorial_map_copy_test.cpp ${hfiles})
target_compile_definitions(Combinatorial_map_copy_test_index PUBLIC USE_COMPACT_CONTAINER_WITH_INDEX)
target_link_libraries(Combinatorial_map_copy_test_index PUBLIC CGAL CGAL::Data)
cgal_add_compilation_test(Combinatorial_map_copy_test_index)

create_single_source_cgal_program(cmap_test_split_attribute.cpp)

# Link with OpenMesh if possible
find_package(OpenMesh QUIET)
if(TARGET OpenMesh::OpenMesh)
  message(STATUS "Found OpenMesh")

  target_link_libraries(Combinatorial_map_copy_test PRIVATE OpenMesh::OpenMesh)
  target_compile_definitions(Combinatorial_map_copy_test PRIVATE -DCGAL_USE_OPENMESH)
  target_link_libraries(Combinatorial_map_copy_test_index PRIVATE OpenMesh::OpenMesh)
  target_compile_definitions(Combinatorial_map_copy_test_index PRIVATE -DCGAL_USE_OPENMESH)
else()
  message(STATUS "NOTICE: Tests will not use OpenMesh.")
endif()
