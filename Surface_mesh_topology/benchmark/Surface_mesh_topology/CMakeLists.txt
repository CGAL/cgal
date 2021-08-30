project(Surface_mesh_topology_Benchmarks)

cmake_minimum_required(VERSION 3.1...3.20)

find_package(CGAL)

# add_definitions(-DCGAL_TRACE_PATH_TESTS)
# add_definitions(-DCGAL_TRACE_CMAP_TOOLS)

set(SOURCE_FILES quadrangulation_computation_benchmarks.cpp
                 path_homotopy_with_schema.cpp)

foreach(cppfile ${SOURCE_FILES})
  create_single_source_cgal_program("${cppfile}")
endforeach()
