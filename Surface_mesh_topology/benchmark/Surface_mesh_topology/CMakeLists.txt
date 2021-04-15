project(Surface_mesh_topology_Benchmarks)

cmake_minimum_required(VERSION 3.1)

find_package(CGAL)

# add_definitions(-DCGAL_TRACE_PATH_TESTS)
# add_definitions(-DCGAL_TRACE_CMAP_TOOLS)

if ( NOT CGAL_FOUND )

  message(FATAL_ERROR "This program requires the CGAL library, and will not be compiled.")

endif()

set(SOURCE_FILES
  quadrangulation_computation_benchmarks.cpp
  path_homotopy_with_schema.cpp
  )

foreach(cppfile ${SOURCE_FILES})
  create_single_source_cgal_program("${cppfile}")
endforeach()
